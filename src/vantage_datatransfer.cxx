#include "../include/vantage_datatransfer.hxx"
#include "bout/bout.hxx"
#include <bout/assert.hxx>
#include <bout/field2d.hxx>
#include <cstddef>
#include <vector>

using namespace NESO::Particles;

// VANTAGE data transfer implementation
// ------------------------------------------------------------------------------
VantageDataTransfer::VantageDataTransfer(
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& project_eval_dg0,
    std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0>& mesh_coupler, Mesh* bout_mesh,
    size_t ndim_vector)
    : neso_mesh(neso_mesh), project_eval_dg0(project_eval_dg0),
      mesh_coupler(mesh_coupler), bout_mesh(bout_mesh),
      dof_kinetic_mesh_scalar(
          std::vector<REAL>(static_cast<size_t>(neso_mesh->get_cell_count()))),
      ndim_vector(ndim_vector),
      dof_kinetic_mesh_vector(std::vector<REAL>(
          ndim_vector * static_cast<size_t>(neso_mesh->get_cell_count()))) {
  // demand that mesh_coupler is not a null pointer
  ASSERT1(mesh_coupler != nullptr)
  // local number of BOUT++ x cells, excluding guards
  const int Nx = bout_mesh->xend - bout_mesh->xstart + 1;
  // local number of BOUT++ y cells, excluding guards
  const int Ny = bout_mesh->yend - bout_mesh->ystart + 1;
  // Get the number of cells in the bout (plasma) mesh owned on this process, excluding guard cells
  num_cells_owned_bout_mesh = static_cast<size_t>(Nx * Ny);
  // a vector used to receive scalar BOUT++ data from the kinetic mesh
  dof_bout_mesh_scalar = std::vector<REAL>(num_cells_owned_bout_mesh);
  dof_bout_mesh_vector = std::vector<REAL>(num_cells_owned_bout_mesh * ndim_vector);
}

void VantageDataTransfer::transfer_scalar_to_plasma_mesh(
    std::vector<REAL>& scalar_kinetic_mesh, Field2D& scalar_plasma_mesh) {
  // some ASSERT required here to check bout_mesh the same
  ASSERT1(scalar_kinetic_mesh.size()
          == static_cast<size_t>(this->neso_mesh->get_cell_count()));
  // we need to port data from the kinetic mesh dofs to the dofs expected by BOUT++ in the loop below
  this->mesh_coupler->backward_transfer(scalar_kinetic_mesh, 1, dof_bout_mesh_scalar);
  std::size_t ic = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      scalar_plasma_mesh(ix, iy) = this->dof_bout_mesh_scalar.at(ic);
      ic++;
    }
  }
  // this fills internal guards
  this->bout_mesh->communicate(scalar_plasma_mesh);
  // apply boundary conditions to fill external guards
  // scalar_field_plasma_mesh.applyBoundary();
  // extrapolate -> Neumann
}

void VantageDataTransfer::transfer_scalar_to_kinetic_mesh(
    Field2D& scalar_plasma_mesh, std::vector<REAL>& scalar_kinetic_mesh) {
  // get scalar from plasma mesh into the dummy vector
  std::size_t ic = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      this->dof_bout_mesh_scalar.at(ic) = scalar_plasma_mesh(ix, iy);
      ic++;
    }
  }
  // transfer data from the dummy vector into the kinetic mesh dofs
  ASSERT1(scalar_kinetic_mesh.size()
          == static_cast<size_t>(this->neso_mesh->get_cell_count()));
  // we need to port data from the kinetic mesh dofs to the dofs expected by BOUT++ in the loop below
  this->mesh_coupler->forward_transfer(dof_bout_mesh_scalar, 1, scalar_kinetic_mesh);
}

void VantageDataTransfer::transfer_scalar_to_particle_property(
    std::vector<REAL>& scalar_kinetic_mesh,
    std::shared_ptr<ParticleGroup>& A_particle_group, std::string particle_property) {
  // check dimensions
  ASSERT1(scalar_kinetic_mesh.size()
          == static_cast<size_t>(this->neso_mesh->get_cell_count()))
  // check the particle property exists with the expected number of dimensions
  ASSERT1(A_particle_group->contains_dat(Sym<REAL>(particle_property), 1));
  // set the kinetic mesh property to NESO-Particles internal variables
  this->project_eval_dg0->set_dofs(1, scalar_kinetic_mesh);
  // set the data from internal variables into the weights
  this->project_eval_dg0->evaluate(A_particle_group, Sym<REAL>(particle_property));
}

void VantageDataTransfer::transfer_scalar_to_particle_property(
    Field2D& scalar_plasma_mesh, std::shared_ptr<ParticleGroup>& A_particle_group,
    std::string particle_property) {
  // check dimensions
  // need some ASSERT to check bout_mesh is the same variables
  this->transfer_scalar_to_kinetic_mesh(scalar_plasma_mesh,
                                        this->dof_kinetic_mesh_scalar);
  this->transfer_scalar_to_particle_property(this->dof_kinetic_mesh_scalar,
                                             A_particle_group, particle_property);
}

void VantageDataTransfer::transfer_particle_property_to_scalar(
    std::shared_ptr<ParticleGroup>& A_particle_group, std::string particle_property,
    std::vector<REAL>& scalar_kinetic_mesh) {
  // check dimensions
  ASSERT1(scalar_kinetic_mesh.size()
          == static_cast<size_t>(this->neso_mesh->get_cell_count()))
  // check the particle property exists with the expected number of dimensions
  ASSERT1(A_particle_group->contains_dat(Sym<REAL>(particle_property), 1));
  // set the particle property to NESO-Particles internal variables
  // some ASSERT to check particle property corresponds to a scalar?
  this->project_eval_dg0->project(A_particle_group, Sym<REAL>(particle_property));
  // project to the kinetic dof vector
  project_eval_dg0->get_dofs(1, scalar_kinetic_mesh);
}

void VantageDataTransfer::transfer_particle_property_to_vector(
    std::shared_ptr<ParticleGroup>& A_particle_group, std::string particle_property,
    std::vector<REAL>& vector_kinetic_mesh) {
  // check dimensions
  ASSERT1(vector_kinetic_mesh.size()
          == this->ndim_vector * static_cast<size_t>(this->neso_mesh->get_cell_count()))
  // check the particle property exists with the expected number of dimensions
  ASSERT1(A_particle_group->contains_dat(Sym<REAL>(particle_property),
                                         static_cast<int>(this->ndim_vector)));
  // set the particle property to NESO-Particles internal variables
  // some ASSERT to check particle property corresponds to a vector?
  this->project_eval_dg0->project(A_particle_group, Sym<REAL>(particle_property));
  // project to the kinetic dof vector
  project_eval_dg0->get_dofs(static_cast<int>(ndim_vector), vector_kinetic_mesh);
}

void VantageDataTransfer::transfer_vector_to_particle_property(
    std::vector<REAL>& vector_kinetic_mesh,
    std::shared_ptr<ParticleGroup>& A_particle_group, std::string particle_property) {
  // check dimensions
  ASSERT1(vector_kinetic_mesh.size()
          == this->ndim_vector * static_cast<size_t>(this->neso_mesh->get_cell_count()))
  // check the particle property exists with the expected number of dimensions
  ASSERT1(A_particle_group->contains_dat(Sym<REAL>(particle_property),
                                         static_cast<int>(this->ndim_vector)));
  // set the kinetic mesh property to NESO-Particles internal variables
  this->project_eval_dg0->set_dofs(static_cast<int>(this->ndim_vector),
                                   vector_kinetic_mesh);
  // set the data from internal variables into the weights
  this->project_eval_dg0->evaluate(A_particle_group, Sym<REAL>(particle_property));
}

void VantageDataTransfer::transfer_vector_to_kinetic_mesh(
    std::vector<Field2D>& vector_plasma_mesh, std::vector<REAL>& vector_kinetic_mesh) {
  // check dimensions
  ASSERT1(vector_plasma_mesh.size() == this->ndim_vector);
  std::size_t ic = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      for (size_t idim = 0; idim < this->ndim_vector; idim++) {
        const size_t jc = ic * this->ndim_vector + idim;
        this->dof_bout_mesh_vector.at(jc) = (vector_plasma_mesh.at(idim))(ix, iy);
      }
      ic++;
    }
  }
  // transfer data from the dummy vector into the kinetic mesh dofs
  ASSERT1(vector_kinetic_mesh.size()
          == ndim_vector * static_cast<size_t>(this->neso_mesh->get_cell_count()));
  // we need to port data from the kinetic mesh dofs to the dofs expected by BOUT++ in the loop below
  this->mesh_coupler->forward_transfer(this->dof_bout_mesh_vector,
                                       static_cast<int>(this->ndim_vector),
                                       vector_kinetic_mesh);
}
void VantageDataTransfer::transfer_vector_to_particle_property(
    std::vector<Field2D>& vector_plasma_mesh,
    std::shared_ptr<ParticleGroup>& A_particle_group, std::string particle_property) {
  this->transfer_vector_to_kinetic_mesh(vector_plasma_mesh,
                                        this->dof_kinetic_mesh_vector);
  this->transfer_vector_to_particle_property(this->dof_kinetic_mesh_vector,
                                             A_particle_group, particle_property);
}

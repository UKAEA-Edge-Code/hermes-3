#pragma once
#include "bout/bout.hxx"
#include <bout/bout_types.hxx>
#include <neso_particles.hpp>
#include <vector>

using namespace NESO::Particles;

/// @brief  Class to handle transfer of data between kinetic mesh and plasma grid in VANTAGE.
class VantageDataTransfer {
public:
  VantageDataTransfer(std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& project_eval_dg0,
    std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0>& mesh_coupler,
    std::shared_ptr<ParticleGroup>& A_particle_group,
    Mesh* bout_mesh, size_t ndim_vector);

  void transfer_scalar_to_plasma_grid(std::vector<REAL>& scalar_kinetic_mesh, Field2D& scalar_plasma_grid);
  void transfer_scalar_to_kinetic_mesh(Field2D& scalar_plasma_grid, std::vector<REAL>& scalar_kinetic_mesh);
  void transfer_scalar_to_particle_property(std::vector<REAL>& scalar_kinetic_mesh, std::string particle_property);
  void transfer_scalar_to_particle_property(Field2D& scalar_plasma_grid, std::string particle_property);
  void transfer_particle_property_to_scalar(std::string particle_property, std::vector<REAL>& scalar_kinetic_mesh);
  void transfer_particle_property_to_vector(std::string particle_property, std::vector<REAL>& vector_kinetic_mesh);

private:
  // internal variables needed for data transfer
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh;
  std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0;
  std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0> mesh_coupler;
  std::shared_ptr<ParticleGroup> A_particle_group;

  // the bout_mesh variable needed for Hermes-3/BOUT++ diagnostics
  Mesh* bout_mesh;

  // variable used to transfer dofs from kinetic
  // to BOUT++ meshes
  size_t num_cells_owned_bout_grid;
  std::vector<REAL> dof_bout_grid_scalar;
  std::vector<REAL> dof_kinetic_mesh_scalar;
  // number of physics vector components
  size_t ndim_vector;
  std::vector<REAL> dof_bout_grid_vector;
  std::vector<REAL> dof_kinetic_mesh_vector;
};

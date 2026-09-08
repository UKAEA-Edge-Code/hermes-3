#include "bout/bout.hxx"
#include "../include/vantage_diagnostics.hxx"

using namespace NESO::Particles;

// VANTAGE diagnostics manager implementation
// ------------------------------------------------------------------------------
VantageDiagnosticsManager::VantageDiagnosticsManager(
    std::string vtkhdf_filename,
    std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0,
    std::shared_ptr<ParticleGroup> A_particle_group,
    BoutReal N_w, BoutReal mass)
    : vtkhdf_filename(vtkhdf_filename),
    neso_mesh(neso_mesh),
    project_eval_dg0(project_eval_dg0),
    A_particle_group(A_particle_group),
    N_w(N_w),
    mass(mass) {}

// Functions for diagnostics on the kinetic mesh
void VantageDiagnosticsManager::write_kinetic_velocity_moment_diagnostics(){
  // get the necessary inputs from the class
  std::string vtkhdf_filename = this->vtkhdf_filename;
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh = this->neso_mesh;
  std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0 = this->project_eval_dg0;
  std::shared_ptr<ParticleGroup> A_particle_group = this->A_particle_group;
  BoutReal N_w = this->N_w;
  BoutReal mass = this->mass;
  // update the necessary particle properties for the moments
  // define the lambda updating the moments
  const size_t ndimv = 2; // number of velocity dimensions
  auto lambda_update_moment_kernels =
        [=](ParticleSubGroupSharedPtr aa) -> void {
      particle_loop(
          "update_moment_kernels", aa,
          [=](auto VELOCITY, auto WEIGHT, auto WEIGHT_V2, auto WEIGHT_V) {
              WEIGHT_V2.at(0) = WEIGHT.at(0) * (VELOCITY.at(0) * VELOCITY.at(0) + VELOCITY.at(1) * VELOCITY.at(1));
              for (int dim = 0; dim < static_cast<int>(ndimv); dim++){
                WEIGHT_V.at(dim) = WEIGHT.at(0) * VELOCITY.at(dim);
              }
          },
          Access::read(Sym<REAL>("VELOCITY")),
          Access::read(Sym<REAL>("WEIGHT")),
          Access::write(Sym<REAL>("WEIGHT_V2")),
          Access::write(Sym<REAL>("WEIGHT_V")))
          ->execute();
    };
  // call the particle loop
  lambda_update_moment_kernels(static_particle_sub_group(A_particle_group));
  // write the data
  VTK::VTKHDF vtk_writer(vtkhdf_filename, neso_mesh->get_comm());
  const std::size_t num_cells_owned_kinetic_mesh = static_cast<size_t>(neso_mesh->get_cell_count());
  // vectors to hold the diagnosed moments
  std::vector<REAL> density(num_cells_owned_kinetic_mesh);
  std::vector<REAL> energy(num_cells_owned_kinetic_mesh);
  std::vector<REAL> gamma(ndimv*num_cells_owned_kinetic_mesh);
  std::vector<REAL> uvector(ndimv*num_cells_owned_kinetic_mesh);
  std::vector<REAL> pressure(num_cells_owned_kinetic_mesh);
  std::vector<REAL> temperature(num_cells_owned_kinetic_mesh);
  // mesh data only CellData not yet filled on each cell
  std::vector<VTK::UnstructuredCell> dvtk0 = neso_mesh->dmh->get_vtk_cell_data();
  std::vector<std::map<std::string, double>> cell_data(num_cells_owned_kinetic_mesh);
  // attempt to explore diagnostics on the DMPlex
  // extract density
  // get data from averages over the diagnostic particle weights
  project_eval_dg0->project(A_particle_group, Sym<REAL>("WEIGHT"));
  // project to the kinetic dof vector
  project_eval_dg0->get_dofs(1, density);
  // energy
  project_eval_dg0->project(A_particle_group, Sym<REAL>("WEIGHT_V2"));
  // project to the kinetic dof vector
  project_eval_dg0->get_dofs(1, energy);
  // mean flow Gamma = nu
  project_eval_dg0->project(A_particle_group, Sym<REAL>("WEIGHT_V"));
  // project to the kinetic dof vector
  project_eval_dg0->get_dofs(2, gamma);
  // scalar variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    // multiply by any factors not handled in the project step
    density.at(ic) *= N_w;
    // (weight factor * mass / 2)
    energy.at(ic) *= 0.5*N_w*mass;
    // insert map entries at this ic
    cell_data.at(ic).insert({"density", density.at(ic)});
    cell_data.at(ic).insert({"energy", energy.at(ic)});
  }
  // vector variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    for (size_t dim=0; dim < ndimv; dim++){
      const size_t jc = ic*ndimv + dim; // compound index covering all cells and dimensions
      // multiply by any factors not handled in the project step
      gamma.at(jc) *= N_w;
      // obtain the derived quantity uvector
      uvector.at(jc) = gamma.at(jc) / density.at(ic);
      // insert a map entry at this ic
      cell_data.at(ic).insert({fmt::format("gamma_{}",dim), gamma.at(jc)});
      cell_data.at(ic).insert({fmt::format("uvector_{}",dim), uvector.at(jc)});
    }
  }
  // derived scalar variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    // calculate pressure = (2E - m n u^2 ) / ndimv
    pressure.at(ic) = 2.0*energy.at(ic);
    for (size_t dim=0; dim < ndimv; dim++){
      const size_t jc = ic*ndimv + dim;
      pressure.at(ic) -= mass*density.at(ic)*uvector.at(jc)*uvector.at(jc);
    }
    pressure.at(ic) /= ndimv;
    temperature.at(ic) = pressure.at(ic)/density.at(ic);
    // insert map entries at this ic
    cell_data.at(ic).insert({"pressure", pressure.at(ic)});
    cell_data.at(ic).insert({"temperature", temperature.at(ic)});
  }
  for (size_t ic=0; ic < static_cast<size_t>(num_cells_owned_kinetic_mesh); ic++){
    // fill the VTK::UnstructuredCell value appropriately
    dvtk0.at(ic).cell_data = cell_data.at(ic);
  }
  vtk_writer.write(dvtk0);
  vtk_writer.close();
}


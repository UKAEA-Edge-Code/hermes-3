#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include <bout/assert.hxx>
#include <bout/field2d.hxx>
#include "../include/component.hxx"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include "../include/vantage_datatransfer.hxx"
#include "../include/vantage_diagnostics.hxx"

using namespace NESO::Particles;

// helper functions for diagnostics

REAL calculate_total_mass(std::vector<REAL>& density,
                     std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh) {
  ASSERT1(density.size() == static_cast<size_t>(neso_mesh->get_cell_count()));
  REAL local_mass = 0.0;
  REAL total_mass = 0.0;
  for (size_t ic = 0;ic < static_cast<size_t>(neso_mesh->get_cell_count()); ic++) {
    local_mass += density.at(ic) * neso_mesh->dmh->get_cell_volume(static_cast<int>(ic));
  }
  MPICHK(
      MPI_Allreduce(&local_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  return total_mass;
}

// helper function to initialise the plasma grid (BOUT++ mesh) diagnostics
Options initialise_plasma_grid_diagnostics(Options& units, Mesh* bout_mesh,
                      //  Field2D& neutral_density,
                      //  Field2D& ion_density,
                      std::string vantage_dump_filepath) {
  // Options object to use to write out diagnostic data of fluid quantities

  const BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  const BoutReal Tnorm = get<BoutReal>(units["eV"]);
  const BoutReal Omega_ci = 1 / get<BoutReal>(units["seconds"]);
  const BoutReal rho_s0 = get<BoutReal>(units["meters"]);
  const BoutReal Bnorm = get<BoutReal>(units["Tesla"]);
  const BoutReal Cs0 = get<BoutReal>(units["meters"])
             / get<BoutReal>(units["seconds"]);

  Options bout_output_data;
  // Add metadata from mesh, e.g. branch cuts
  bout_mesh->outputVars(bout_output_data);
  // Add Rxy, Zxy coordinate data
  Field2D Rxy;
  Field2D Rxy_corners;
  Field2D Rxy_lower_right_corners;
  Field2D Rxy_upper_right_corners;
  Field2D Rxy_upper_left_corners;
  Field2D Zxy;
  Field2D Zxy_corners;
  Field2D Zxy_lower_right_corners;
  Field2D Zxy_upper_right_corners;
  Field2D Zxy_upper_left_corners;
  // mesh->get(ivertex, "ivertex_lower_left_corners");
  bout_mesh->get(Rxy, "Rxy");
  bout_mesh->get(Rxy_corners, "Rxy_corners");
  bout_mesh->get(Rxy_lower_right_corners, "Rxy_lower_right_corners");
  bout_mesh->get(Rxy_upper_right_corners, "Rxy_upper_right_corners");
  bout_mesh->get(Rxy_upper_left_corners, "Rxy_upper_left_corners");
  bout_mesh->get(Zxy, "Zxy");
  bout_mesh->get(Zxy_corners, "Zxy_corners");
  bout_mesh->get(Zxy_lower_right_corners, "Zxy_lower_right_corners");
  bout_mesh->get(Zxy_upper_right_corners, "Zxy_upper_right_corners");
  bout_mesh->get(Zxy_upper_left_corners, "Zxy_upper_left_corners");
  set_with_attrs(bout_output_data["Rxy"], Rxy, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Rxy_corners"], Rxy_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Rxy_lower_right_corners"], Rxy_lower_right_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Rxy_upper_right_corners"], Rxy_upper_right_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Rxy_upper_left_corners"], Rxy_upper_left_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Zxy"], Zxy, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Zxy_corners"], Zxy_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Zxy_lower_right_corners"], Zxy_lower_right_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Zxy_upper_right_corners"], Zxy_upper_right_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["Zxy_upper_left_corners"], Zxy_upper_left_corners, {
      {"units", "m"},
      {"conversion", 1}, // Already in SI units
    });
  set_with_attrs(bout_output_data["y_boundary_guards"], 2, {
      {"source", "vantage -- should be provided by BOUT++"}
    });

  // Add metadata with normalisation factors
  set_with_attrs(bout_output_data["Tnorm"], Tnorm,
                 {{"units", "eV"},
                  {"conversion", 1}, // Already in SI units
                  {"standard_name", "temperature normalisation"},
                  {"long_name", "temperature normalisation"}});
  set_with_attrs(bout_output_data["Nnorm"], Nnorm,
                 {{"units", "m^-3"},
                  {"conversion", 1},
                  {"standard_name", "density normalisation"},
                  {"long_name", "Number density normalisation"}});
  set_with_attrs(bout_output_data["Bnorm"], Bnorm,
                 {{"units", "T"},
                  {"conversion", 1},
                  {"standard_name", "magnetic field normalisation"},
                  {"long_name", "Magnetic field normalisation"}});
  set_with_attrs(bout_output_data["Cs0"], Cs0,
                 {{"units", "m/s"},
                  {"conversion", 1},
                  {"standard_name", "velocity normalisation"},
                  {"long_name", "Sound speed normalisation"}});
  set_with_attrs(bout_output_data["Omega_ci"], Omega_ci,
                 {{"units", "s^-1"},
                  {"conversion", 1},
                  {"standard_name", "frequency normalisation"},
                  {"long_name", "Cyclotron frequency normalisation"}});
  set_with_attrs(bout_output_data["rho_s0"], rho_s0,
                 {{"units", "m"},
                  {"conversion", 1},
                  {"standard_name", "length normalisation"},
                  {"long_name", "Gyro-radius length normalisation"}});

  bout::OptionsIO::create(vantage_dump_filepath)->write(bout_output_data);
  return bout_output_data;
}

// VANTAGE diagnostics manager implementation
// ------------------------------------------------------------------------------
VantageDiagnosticsManager::VantageDiagnosticsManager(
    std::string vtkhdf_filename,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    std::shared_ptr<ParticleGroup>& A_particle_group,
    std::shared_ptr<VantageDataTransfer>& data_transfer,
    BoutReal N_w, BoutReal mass, Mesh* bout_mesh,
    Options& units, std::string vantage_dump_filepath)
    : vtkhdf_filename(vtkhdf_filename),
    neso_mesh(neso_mesh),
    A_particle_group(A_particle_group),
    data_transfer(data_transfer),
    N_w(N_w),
    mass(mass),
    bout_mesh(bout_mesh),
    units(units),
    vantage_dump_filepath(vantage_dump_filepath) {
      // initialise vectors for storing the moments on the kinetic mesh
      const size_t ndimv = this->ndimv;
      const std::size_t num_cells_owned_kinetic_mesh = static_cast<size_t>(neso_mesh->get_cell_count());
      density = std::vector<REAL>(num_cells_owned_kinetic_mesh);
      energy = std::vector<REAL>(num_cells_owned_kinetic_mesh);
      gamma = std::vector<REAL>(ndimv*num_cells_owned_kinetic_mesh);
      uvector = std::vector<REAL>(ndimv*num_cells_owned_kinetic_mesh);
      pressure = std::vector<REAL>(num_cells_owned_kinetic_mesh);
      temperature = std::vector<REAL>(num_cells_owned_kinetic_mesh);
      // initialise BOUT++ diagnostic on BOUT++ mesh
      // Object for VANTAGE dump files
      vantage_dump_writer =
          bout::OptionsIO::create({{"file", vantage_dump_filepath}, {"append", true}});
      bout_output_data =
        initialise_plasma_grid_diagnostics(units, bout_mesh, vantage_dump_filepath);
      // initialise plasma grid variables
      density_plasma_grid = Field2D(0.0, bout_mesh);
      energy_plasma_grid = Field2D(0.0, bout_mesh);
      pressure_plasma_grid = Field2D(0.0, bout_mesh);
      temperature_plasma_grid = Field2D(0.0, bout_mesh);
      // local number of BOUT++ x cells, excluding guards
      const int Nx = bout_mesh->xend - bout_mesh->xstart + 1;
      // local number of BOUT++ y cells, excluding guards
      const int Ny = bout_mesh->yend - bout_mesh->ystart + 1;
      // Get the number of cells in the bout (plasma) mesh owned on this process, excluding guard cells
      const size_t num_cells_owned_bout_mesh = static_cast<size_t>(Nx*Ny);
      // a vector used to receive scalar BOUT++ data from the kinetic mesh
      dof_bout_mesh_scalar = std::vector<REAL>(num_cells_owned_bout_mesh);
    }

// Functions for diagnostics on the kinetic mesh
void VantageDiagnosticsManager::update_kinetic_velocity_moments(){
  // get the necessary inputs from the class
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh = this->neso_mesh;
  std::shared_ptr<ParticleGroup> A_particle_group = this->A_particle_group;
  BoutReal N_w = this->N_w;
  BoutReal mass = this->mass;
  // get the necessary private diagnostic variables
  // vectors to hold the diagnosed moments
  // use references here since we want to update the members
  std::vector<REAL>& density = this->density;
  std::vector<REAL>& energy = this->energy;
  std::vector<REAL>& gamma = this->gamma;
  std::vector<REAL>& uvector = this->uvector;
  std::vector<REAL>& pressure = this->pressure;
  std::vector<REAL>& temperature = this->temperature;

  // update the necessary particle properties for the moments
  // define the lambda updating the moments
  const size_t ndimv = this->ndimv; // number of velocity dimensions
  const std::size_t num_cells_owned_kinetic_mesh = static_cast<size_t>(neso_mesh->get_cell_count());
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
  // extract density
  this->data_transfer->transfer_particle_property_to_scalar("WEIGHT", density);
  // energy
  this->data_transfer->transfer_particle_property_to_scalar("WEIGHT_V2", energy);
  // mean flow Gamma = nu
  this->data_transfer->transfer_particle_property_to_vector("WEIGHT_V", gamma);
  // scalar variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    // multiply by any factors not handled in the project step
    density.at(ic) *= N_w;
    // (weight factor * mass / 2)
    energy.at(ic) *= 0.5*N_w*mass;
  }
  // vector variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    for (size_t dim=0; dim < ndimv; dim++){
      const size_t jc = ic*ndimv + dim; // compound index covering all cells and dimensions
      // multiply by any factors not handled in the project step
      gamma.at(jc) *= N_w;
      // obtain the derived quantity uvector
      uvector.at(jc) = gamma.at(jc) / density.at(ic);
    }
  }
  // derived scalar variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    // calculate pressure = (2E - m n u^2 ) / ndimv
    pressure.at(ic) = 2.0*energy.at(ic);
    for (size_t dim=0; dim < ndimv; dim++){
      const size_t jc = (ic*ndimv) + dim;
      pressure.at(ic) -= mass*density.at(ic)*uvector.at(jc)*uvector.at(jc);
    }
    pressure.at(ic) /= static_cast<REAL>(ndimv);
    temperature.at(ic) = pressure.at(ic)/density.at(ic);
  }
}

// Function to save a VTKHDF file, writing the private member
// velocity moments and the mesh in VTK compatible format
// we pass in the ion_density here to enable testing of mass conservation
void VantageDiagnosticsManager::write_kinetic_velocity_moment_diagnostics(int istep,
        std::vector<REAL>& ion_density){
  // get the necessary inputs from the class
  const std::string vtkhdf_filename = fmt::format("{}.istep.{}.vtkhdf",this->vtkhdf_filename,istep);
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh = this->neso_mesh;

  // write the data
  VTK::VTKHDF vtk_writer(vtkhdf_filename, neso_mesh->get_comm());
  const std::size_t num_cells_owned_kinetic_mesh = static_cast<size_t>(neso_mesh->get_cell_count());
  // vectors to hold the diagnosed moments
  std::vector<REAL> density = this->density;
  std::vector<REAL> energy = this->energy;
  std::vector<REAL> gamma = this->gamma;
  std::vector<REAL> uvector = this->uvector;
  std::vector<REAL> pressure = this->pressure;
  std::vector<REAL> temperature = this->temperature;
  // mesh data only CellData not yet filled on each cell
  std::vector<VTK::UnstructuredCell> dvtk0 = neso_mesh->dmh->get_vtk_cell_data();
  std::vector<std::map<std::string, double>> cell_data(num_cells_owned_kinetic_mesh);
  // scalar variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    // insert map entries at this ic
    cell_data.at(ic).insert({"density", density.at(ic)});
    cell_data.at(ic).insert({"energy", energy.at(ic)});
    cell_data.at(ic).insert({"pressure", pressure.at(ic)});
    cell_data.at(ic).insert({"temperature", temperature.at(ic)});
    // write cell volume for convenience in later post-processing analysis
    cell_data.at(ic).insert({"cellvolume", neso_mesh->dmh->get_cell_volume(static_cast<int>(ic))});
    // write the "ion density" for testing purposes only
    cell_data.at(ic).insert({"ion_density", ion_density.at(ic)});
  }
  // vector variables
  for (size_t ic=0; ic < num_cells_owned_kinetic_mesh; ic++){
    for (size_t dim=0; dim < ndimv; dim++){
      const size_t jc = ic*ndimv + dim; // compound index covering all cells and dimensions
      // insert a map entry at this ic
      cell_data.at(ic).insert({fmt::format("gamma_{}",dim), gamma.at(jc)});
      cell_data.at(ic).insert({fmt::format("uvector_{}",dim), uvector.at(jc)});
    }
  }
  for (size_t ic=0; ic < static_cast<size_t>(num_cells_owned_kinetic_mesh); ic++){
    // fill the VTK::UnstructuredCell value appropriately
    dvtk0.at(ic).cell_data = cell_data.at(ic);
  }
  vtk_writer.write(dvtk0);
  vtk_writer.close();
}

void VantageDiagnosticsManager::transfer_moments_to_plasma_grid(){
  this->data_transfer->transfer_scalar_to_plasma_grid(
    this->density,
    this->density_plasma_grid);
  this->data_transfer->transfer_scalar_to_plasma_grid(
    this->energy,
    this->energy_plasma_grid);
  this->data_transfer->transfer_scalar_to_plasma_grid(
    this->pressure,
    this->pressure_plasma_grid);
  this->data_transfer->transfer_scalar_to_plasma_grid(
    this->temperature,
    this->temperature_plasma_grid);
}

void VantageDiagnosticsManager::write_bout_diagnostics(std::vector<REAL>& ion_density_kinetic_mesh, Field2D& Siz,
                        Field2D& Srec,
                        BoutReal particle_time) {

  Field2D ion_density = Field2D{0.0, this->bout_mesh};
  this->data_transfer->transfer_scalar_to_plasma_grid(
    ion_density_kinetic_mesh, ion_density);
  // extract the units
  const BoutReal Nnorm = get<BoutReal>(this->units["inv_meters_cubed"]);
  // const BoutReal Tnorm = get<BoutReal>(units["eV"]);
  const BoutReal Omega_ci = 1 / get<BoutReal>(this->units["seconds"]);
  // const BoutReal rho_s0 = get<BoutReal>(units["meters"]);
  // const BoutReal Bnorm = get<BoutReal>(units["Tesla"]);
  // const BoutReal Cs0 = get<BoutReal>(units["meters"])
  //            / get<BoutReal>(units["seconds"]);
  Field2D neutral_density = this->density_plasma_grid;
  set_with_attrs(this->bout_output_data["neutral_density"], neutral_density,
                 {{"time_dimension", "t"}});

  set_with_attrs(this->bout_output_data["ion_density"], ion_density, {{"time_dimension", "t"}});

  set_with_attrs(this->bout_output_data["Nn"], neutral_density,
                 {{"time_dimension", "t"},
                  {"units", "m^-3"},
                  {"conversion", Nnorm},
                  {"standard_name", "Density"},
                  {"long_name", "Kinetic neutral density"},
                  {"species", "kinetic neutrals"},
                  {"source", "vantage"}});

  set_with_attrs(this->bout_output_data["Siz"], Siz,
                 {{"time_dimension", "t"},
                  {"units", "m^-3 s^-1"},
                  {"conversion", Nnorm * Omega_ci},
                  {"standard_name", "Density source"},
                  {"long_name", "Ionisation density source"},
                  {"species", "kinetic neutrals"},
                  {"source", "vantage"}});

  set_with_attrs(this->bout_output_data["Srec"], Srec,
                 {{"time_dimension", "t"},
                  {"units", "m^-3 s^-1"},
                  {"conversion", Nnorm * Omega_ci},
                  {"standard_name", "Density source"},
                  {"long_name", "Recombination density source"},
                  {"species", "kinetic neutrals"},
                  {"source", "vantage"}});

  set_with_attrs(this->bout_output_data["t_array"], particle_time, {{"time_dimension", "t"}});

  // Append data to file
  this->vantage_dump_writer->write(this->bout_output_data);
  // Ensure buffer is written to disk to avoid crash data loss
  this->vantage_dump_writer->flush();
}

std::vector<REAL> VantageDiagnosticsManager::get_density_kinetic_mesh(){
  return this->density;
}



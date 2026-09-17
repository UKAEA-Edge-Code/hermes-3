#pragma once
#include "bout/bout.hxx"
#include <bout/bout_types.hxx>
#include <memory>
#include <neso_particles.hpp>
#include <string>
#include <vector>
#include "../include/vantage_datatransfer.hxx"
#include "../include/vantage_sources.hxx"

using namespace NESO::Particles;

REAL calculate_total_mass(std::vector<REAL>& density,
          std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh);
REAL calculate_total_mass(Field2D& density,
          std::vector<REAL>& neso_cell_volume_on_bout_mesh);

/// @brief  Class to manage diagnostics from VANTAGE.
class VantageDiagnosticsManager {
public:
  VantageDiagnosticsManager(std::string vtkhdf_filename,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    // neso_mesh cell volumes on the BOUT++ mesh
    std::vector<REAL>& neso_cell_volumes,
    std::shared_ptr<ParticleGroup>& A_particle_group,
    std::shared_ptr<VantageDataTransfer>& data_transfer,
    std::shared_ptr<VantageSourceManager>& source_manager,
    BoutReal N_w, BoutReal mass,
    Mesh* bout_mesh, Options& units,
    std::string vantage_dump_filepath);

  // compute the kinetic velocity moments and
  // store in private variables
  void update_kinetic_velocity_moments();
  // write kinetic diagnostics to a vtkhdf file
  void write_kinetic_velocity_moment_diagnostics(int istep, std::vector<REAL>& ion_density);
  // transfer kinetic moments to BOUT++ grid
  void transfer_moments_to_plasma_grid();
  // write BOUT++ style diagnostics on the BOUT++ grid
  void write_bout_diagnostics(
          Field2D& ion_density,
          BoutReal particle_time);
  std::vector<REAL> get_density_kinetic_mesh();
  // Field2D transfer_scalar_to_plasma_grid(std::vector<REAL>& scalar_field);

private:
  // internal variables needed for diagnostics
  std::string vtkhdf_filename;
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh;
  std::shared_ptr<ParticleGroup> A_particle_group;
  std::shared_ptr<VantageDataTransfer> data_transfer;
  // pointer to vantage source_manager for diagnostics
  std::shared_ptr<VantageSourceManager> source_manager;
  BoutReal N_w;
  BoutReal mass;

  // the bout_mesh variable needed for Hermes-3/BOUT++ diagnostics
  Mesh* bout_mesh;
  Options& units;
  std::string vantage_dump_filepath;
  Options bout_output_data; // Options object to hold output data for VANTAGE diagnostics
  std::unique_ptr<bout::OptionsIO>
      vantage_dump_writer; // OptionsIO object to write VANTAGE diagnostics

  // variables used to store the moments of
  // the neutral distribution function, on
  // the kinetic mesh
  const size_t ndimv = 2; // number of velocity dimensions
  std::vector<REAL> density;
  std::vector<REAL> energy;
  std::vector<REAL> gamma;
  std::vector<REAL> uvector;
  std::vector<REAL> pressure;
  std::vector<REAL> temperature;

  // variable used to transfer dofs from kinetic
  // to BOUT++ meshes
  std::vector<REAL> dof_bout_mesh_scalar;
  // variables used to store the moments of the
  // neutral distribution function projected on
  // to the BOUT++ grid
  Field2D density_plasma_grid;
  Field2D energy_plasma_grid;
  Field2D pressure_plasma_grid;
  Field2D temperature_plasma_grid;
  // n.b. only treat scalar variables for now
};

#pragma once
#include "bout/bout.hxx"
#include <bout/bout_types.hxx>
#include <neso_particles.hpp>
#include <string>
#include <vector>

using namespace NESO::Particles;

REAL calculate_total_mass(std::vector<REAL>& density,
          std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh);

/// @brief  Class to manage diagnostics from VANTAGE.
class VantageDiagnosticsManager {
public:
  VantageDiagnosticsManager(std::string vtkhdf_filename,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& project_eval_dg0,
    std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0>& mesh_coupler,
    std::shared_ptr<ParticleGroup>& A_particle_group,
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
          std::vector<REAL>& ion_density_kinetic_mesh,
          Field2D& Siz,
          Field2D& Srec,
          BoutReal particle_time);
  std::vector<REAL> get_density_kinetic_mesh();
  Field2D transfer_scalar_to_plasma_grid(std::vector<REAL>& scalar_field);

private:
  // internal variables needed for diagnostics
  std::string vtkhdf_filename;
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh;
  std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0;
  std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0> mesh_coupler;
  std::shared_ptr<ParticleGroup> A_particle_group;
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

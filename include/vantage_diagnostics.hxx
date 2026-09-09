#pragma once
#include "bout/bout.hxx"
#include <bout/bout_types.hxx>
#include <neso_particles.hpp>
#include <string>
#include <vector>

using namespace NESO::Particles;

BoutReal calculate_total_mass(Field2D& density,
          std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh);

/// @brief  Class to manage diagnostics from VANTAGE.
class VantageDiagnosticsManager {
public:
  VantageDiagnosticsManager(std::string vtkhdf_filename,
    std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0,
    std::shared_ptr<ParticleGroup> A_particle_group,
    BoutReal N_w, BoutReal mass,
    Mesh* bout_mesh, Options& units,
    std::string vantage_dump_filepath);

  // compute the kinetic velocity moments and
  // store in private variables
  void update_kinetic_velocity_moments();
  // write kinetic diagnostics to a vtkhdf file
  void write_kinetic_velocity_moment_diagnostics();
  // write BOUT++ style diagnostics on the BOUT++ grid
  void write_bout_diagnostics(
          Field2D& neutral_density,
          Field2D& ion_density,
          Field2D& Siz,
          Field2D& Srec,
          BoutReal particle_time);
private:
  // internal variables needed for diagnostics
  std::string vtkhdf_filename;
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh;
  std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0;
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
};

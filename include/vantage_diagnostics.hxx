#pragma once
#include "bout/bout.hxx"
#include <neso_particles.hpp>
#include <vector>

using namespace NESO::Particles;


/// @brief  Class to manage diagnostics from VANTAGE.
class VantageDiagnosticsManager {
public:
  VantageDiagnosticsManager(std::string vtkhdf_filename,
    std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0,
    std::shared_ptr<ParticleGroup> A_particle_group,
    BoutReal N_w, BoutReal mass);

  // compute the kinetic velocity moments and
  // store in private variables
  void update_kinetic_velocity_moments();
  // write kinetic diagnostics to a vtkhdf file
  void write_kinetic_velocity_moment_diagnostics();

private:
  // internal variables needed for diagnostics
  std::string vtkhdf_filename;
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh;
  std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0;
  std::shared_ptr<ParticleGroup> A_particle_group;
  BoutReal N_w;
  BoutReal mass;

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

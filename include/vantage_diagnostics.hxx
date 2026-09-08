#pragma once
#include "bout/bout.hxx"
#include <neso_particles.hpp>

using namespace NESO::Particles;


/// @brief  Class to manage diagnostics from VANTAGE.
class VantageDiagnosticsManager {
public:
  VantageDiagnosticsManager(std::string vtkhdf_filename,
    std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG> project_eval_dg0,
    std::shared_ptr<ParticleGroup> A_particle_group,
    BoutReal N_w, BoutReal mass);

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
};

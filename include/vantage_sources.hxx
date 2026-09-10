#pragma once
#include "bout/bout.hxx"
#include <neso_particles.hpp>
#include <reactions_lib/common_transformations.hpp>
#include <reactions_lib/transformation_wrapper.hpp>

using namespace NESO::Particles;
using namespace VANTAGE::Reactions;

/// @brief Data struct to hold information about a reaction source.
/// @param reaction_name Name of the reaction, e.g. "ionistaion"
/// @param source_name Name of the source, e.g. Siz (ion density source due to
/// ionisation).
/// @param accumulator CellwiseAccumulator to use to accumulate the source term for this
/// reaction.
/// @param particle_group ParticleGroup to which this source applies.
/// @param zeroer TransformationStrategy to use to zero the source term dat after
/// accumulation.
struct VantageSource {
  std::string hermes_source_name;
  std::string vantage_source_name;
  std::shared_ptr<CellwiseAccumulator<REAL>> accumulator;
  std::shared_ptr<ParticleGroup> particle_group;
  std::shared_ptr<TransformationStrategy> zeroer;
  Field2D source_data_plasma_grid;
};

/// @brief  Class to manage reaction channel sources from VANTAGE.
/// Source terms from VANTAGE are extracted from the accumulator.
/// These are then converted to actual sources, e.g. units of m^-3 s^-1 for a
/// density source.
class VantageSourceManager {
public:
  VantageSourceManager(std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                      std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0>& mesh_coupler_dg0,
                      std::vector<REAL>& dof_kinetic_mesh_scalar,
                      std::vector<REAL>& dof_bout_mesh_scalar,
                      Mesh* bout_mesh, Options& units);

  Mesh* bout_mesh;

  // Register new source
  void add_source(const std::string& hermes_source_name,
                  const std::string& vantage_source_name,
                  std::shared_ptr<CellwiseAccumulator<REAL>> accumulator,
                  std::shared_ptr<ParticleGroup> particle_group,
                  std::shared_ptr<TransformationStrategy> zeroer);

  // Update the Hermes-3 source field using the accumulated data from corresponding
  // VANTAGE source
  void update_source(const std::string& hermes_source_name, double dt);

  // Call update_source on all sources
  void update_all_sources(double dt);

  // Return data for a given Hermes-3 source name
  Field2D get_data(const std::string& hermes_source_name);

private:
  std::map<std::string, VantageSource> sources;
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh;
  std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0> mesh_coupler_dg0;
  std::vector<REAL> dof_kinetic_mesh_scalar;
  std::vector<REAL> dof_bout_mesh_scalar;
  Options& units;
};
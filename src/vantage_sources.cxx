#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include <bout/assert.hxx>
#include "../include/component.hxx"
#include <memory>
#include <neso_particles.hpp>
#include "../include/vantage_sources.hxx"
#include <reactions_lib/common_transformations.hpp>
#include <reactions_lib/transformation_wrapper.hpp>
#include <vector>
#include "../include/vantage_datatransfer.hxx"


using namespace NESO::Particles;
using namespace VANTAGE::Reactions;

// VANTAGE source manager implementation
// ------------------------------------------------------------------------------
VantageSourceManager::VantageSourceManager(
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    std::shared_ptr<VantageDataTransfer>& data_transfer,
    Mesh* bout_mesh,
    Options& units)
    : bout_mesh(bout_mesh),
    neso_mesh(neso_mesh),
    data_transfer(data_transfer),
    units(units) {}

// Register new source with the manager and initialise its data
void VantageSourceManager::add_source(
    const std::string& hermes_source_name, const std::string& vantage_source_name,
    std::shared_ptr<CellwiseAccumulator<REAL>> accumulator,
    std::shared_ptr<ParticleGroup> particle_group,
    std::shared_ptr<TransformationStrategy> zeroer) {

  Field2D source_data_plasma_grid{bout_mesh};
  source_data_plasma_grid = 0.0;
  const int num_cells_owned_kinetic_mesh = neso_mesh->get_cell_count();
  std::vector<REAL> source_data_kinetic_mesh(static_cast<size_t>(num_cells_owned_kinetic_mesh));

  VantageSource source{
      hermes_source_name, vantage_source_name,
      accumulator, particle_group, zeroer,
      source_data_plasma_grid,
      source_data_kinetic_mesh};

  this->sources[hermes_source_name] = source;
}

// Return source data on plasma grid
Field2D VantageSourceManager::get_plasma_grid_data(const std::string& hermes_source_name) {
  return this->sources[hermes_source_name].source_data_plasma_grid;
}

// Return source data on kinetic mesh
std::vector<REAL> VantageSourceManager::get_kinetic_mesh_data(const std::string& hermes_source_name) {
  return this->sources[hermes_source_name].source_data_kinetic_mesh;
}

// Update the source from VANTAGE and reset the VANTAGE data/accumulator
void VantageSourceManager::update_source(const std::string& hermes_source_name,
                                         double dt) {

  VantageSource& source = this->sources[hermes_source_name];
  BoutReal N_w = get<BoutReal>(units["N_w"]);

  std::vector<CellData<double>> accumulated_1d =
      source.accumulator->get_cell_data(source.vantage_source_name);
  size_t naccumulated = accumulated_1d.size();
  ASSERT1(naccumulated == source.source_data_kinetic_mesh.size());
  for (size_t ic = 0; ic < naccumulated; ic++){
    source.source_data_kinetic_mesh.at(ic) = accumulated_1d[ic]->at(0, 0)   // Total weight
          * N_w                                                   // Total particles
          / neso_mesh->dmh->get_cell_volume(static_cast<int>(ic)) // Total density
          / dt;                                                   // Density source;
  }
  // copy accumulated data into the BOUT++ Field2D variable
  this->data_transfer->transfer_scalar_to_plasma_grid(
    source.source_data_kinetic_mesh, source.source_data_plasma_grid);
  // Reset the accumulator object
  source.accumulator->zero_buffer(source.vantage_source_name);
  // Reset the accumulated source data on the particle
  source.zeroer->transform(std::make_shared<ParticleSubGroup>(source.particle_group));
}

// Update all sources
void VantageSourceManager::update_all_sources(double dt) {
  for (auto& [hermes_source_name, source] : this->sources) {
    update_source(hermes_source_name, dt);
  }
}

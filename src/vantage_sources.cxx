#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include <bout/assert.hxx>
#include "../include/component.hxx"
#include <neso_particles.hpp>
#include "../include/vantage_sources.hxx"
#include <reactions_lib/common_transformations.hpp>
#include <reactions_lib/transformation_wrapper.hpp>
#include <vector>


using namespace NESO::Particles;
using namespace VANTAGE::Reactions;

// VANTAGE source manager implementation
// ------------------------------------------------------------------------------
VantageSourceManager::VantageSourceManager(
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0>& mesh_coupler_dg0,
    std::vector<REAL>& dof_kinetic_mesh_scalar,
    std::vector<REAL>& dof_bout_mesh_scalar,
    Mesh* bout_mesh,
    Options& units)
    : bout_mesh(bout_mesh),
    neso_mesh(neso_mesh),
    mesh_coupler_dg0(mesh_coupler_dg0),
    dof_kinetic_mesh_scalar(dof_kinetic_mesh_scalar), dof_bout_mesh_scalar(dof_bout_mesh_scalar),
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

// Return source data
Field2D VantageSourceManager::get_data(const std::string& hermes_source_name) {
  return this->sources[hermes_source_name].source_data_plasma_grid;
}

// Update the source from VANTAGE and reset the VANTAGE data/accumulator
void VantageSourceManager::update_source(const std::string& hermes_source_name,
                                         double dt) {

  VantageSource& source = this->sources[hermes_source_name];
  BoutReal N_w = get<BoutReal>(units["N_w"]);

  std::vector<CellData<double>> accumulated_1d =
      source.accumulator->get_cell_data(source.vantage_source_name);
  size_t naccumulated = accumulated_1d.size();
  if (mesh_coupler_dg0 != nullptr){
    ASSERT1(source.source_data_kinetic_mesh.size() > dof_bout_mesh_scalar.size())
    // copy accumulated data into the relevant kinetic dof variable
    for (size_t ic = 0; ic < naccumulated; ic++){
      dof_kinetic_mesh_scalar.at(ic) = accumulated_1d[ic]->at(0, 0);
      source.source_data_kinetic_mesh.at(ic) = dof_kinetic_mesh_scalar.at(ic);
    }
    // use the transform from kinetic to bout mesh
    mesh_coupler_dg0->backward_transfer(dof_kinetic_mesh_scalar, 1, dof_bout_mesh_scalar);
  } else {
    ASSERT1(source.source_data_kinetic_mesh.size() == dof_bout_mesh_scalar.size())
    // copy accumulated data directly into the relevant bout dof variable
    for (size_t ic = 0; ic < naccumulated; ic++){
      dof_bout_mesh_scalar.at(ic) = accumulated_1d[ic]->at(0, 0);
      source.source_data_kinetic_mesh.at(ic) = dof_bout_mesh_scalar.at(ic);
    }
  }
  std::size_t ic = 0;
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      source.source_data_plasma_grid(ix, iy) =
          dof_bout_mesh_scalar.at(ic)                            // Total weight
          * N_w                                                   // Total particles
          / neso_mesh->dmh->get_cell_volume(static_cast<int>(ic)) // Total density
          / dt;                                                   // Density source

      ic++;
    }
  }

  // Fill internal guards
  bout_mesh->communicate(source.source_data_plasma_grid);
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

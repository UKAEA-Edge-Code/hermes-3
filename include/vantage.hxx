#pragma once
#include "../include/component.hxx"
#include "bout/bout.hxx"
#include "bout/petsclib.hxx"
#include <neso_particles.hpp>
#include <neso_rng_toolkit.hpp>
#include <reactions/reactions.hpp>

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
  Field2D source_data;
};

/// @brief  Class to manage reaction channel sources from VANTAGE.
/// Source terms from VANTAGE are extracted from the accumulator.
/// These are then converted to actual sources, e.g. units of m^-3 s^-1 for a
/// density source.
class VantageSourceManager {
public:
  VantageSourceManager(std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
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
  Options& units;
};

// Need to declare empty struct because the Monitor needs it and it must be
// before component construction as it has the monitor as a member.
// See https://bout-dev.readthedocs.io/en/latest/user_docs/time_integration.html#monitoring-the-simulation-output
// for more information on monitors.
struct Vantage;

/// @brief  Monitor to schedule kinetic iterations.
/// @param solver pointer to the solver.
/// @param time normalised sim time provided by BOUT++
/// @param iter current iteration number provided by BOUT++
/// @param nout number of outputs provided by BOUT++
class VantageMonitor : public Monitor {
public:
  explicit VantageMonitor(Vantage* vantage) : vantage(vantage) {}
  int call(Solver* solver, BoutReal time, int iter, int nout) override;

private:
  Vantage* vantage;
};

struct Vantage : public Component {
  Vantage(std::string name, Options& options, Solver* solver);

  ~Vantage(); // Destructor for VANTAGE related cleanup
  void finally(const Options& state) override;
  void transform_impl(GuardedOptions& state) override;
  void outputVars(Options& state) override;

  // Function with the kinetic loop.
  // time is the normalised VANTAGE monitor frequency.
  int advance_vantage(BoutReal time);

private:
  bool test_mass_conservation;
  BoutReal N_w;
  REAL dt;
  int nsteps;
  int num_cells_owned; // Number of VANTAGE cells owned per rank

  Options bout_output_data; // Options object to hold output data for VANTAGE diagnostics
  std::unique_ptr<bout::OptionsIO>
      vantage_dump_writer; // OptionsIO object to write VANTAGE diagnostics

  int mpi_rank;    // Current rank ID
  Mesh* bout_mesh; // Pointer to the BOUT++ mesh object
  Field2D ion_density, neutral_density, total_density;
  Field2D initial_neutral_density; // Initial VANTAGE kinetic neutral density
  BoutReal total_mass_initial, total_mass;
  std::string dmplex_filepath, vantage_dump_filepath,
      particle_data_filepath; // Path for output files

  // Physics
  BoutReal particle_time;
  std::string neutral_species; // Neutral species to simulate with VANTAGE
  std::string ion_species;
  bool plasma_coupling; // Whether to read plasma fields from the state or not
  PetscLib petsc_lib;   // Ensures PETSc is initialized for the lifetime of this component

  // DMPLex/VANTAGE stuff
  DM dm;
  std::shared_ptr<PetscInterface::DMPlexInterface> neso_mesh;
  std::shared_ptr<SYCLTarget> sycl_target;
  std::shared_ptr<PetscInterface::BoundaryInteraction2D>
      b2d; // Boundary interaction object
  std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>
      dg0;                        // DMPlex projection object
  std::shared_ptr<H5Part> h5part; // HDF5 particle output object
  std::vector<REAL>
      h_project1; // Buffer for scalar projection/evaluation of NESO-Particles properties
  std::shared_ptr<ParticleGroup> A_particle_group; // Particle group for main neutrals
  std::shared_ptr<ParticleGroup> marker_group;     // Particle group for rec markers

  // Needed for Vantage::apply_boundary_conditions
  std::shared_ptr<BoundaryReflection> reflection; // Boundary reflection object
  void apply_boundary_conditions(ParticleSubGroupSharedPtr aa);

  // These classes don't have a default constructor so need to be initialised as a unique_ptr
  std::unique_ptr<VantageSourceManager>
      source_manager;           // Manager for VANTAGE reaction sources
  VantageMonitor monitor{this}; // Output monitor to schedule VANTAGE iterations

  std::unique_ptr<ReactionController> reaction_controller;
  std::unique_ptr<ReactionController> recombination_controller;
};

namespace {
RegisterComponent<Vantage> registercomponentvantage("vantage");
}

/**
 * @brief Function to calculate cell volumes.
 *
 * @param sycl_target SYCLTargetSharedPtr to use for communication.
 * @param mesh object.
 *
 */

/**
 * @brief Function to calculate particle positions and velocities from a Maxwellian.
 *
 * @param sycl_target SYCLTargetSharedPtr to use for communication.
 * @param mesh object.
 * @param particle_spec ParticleSpec to use for the returned ParticleSet.
 * @param npart_per_cell Number of particles to create per cell.
 * @param weight Weight to assign to each particle.
 * @param std_dev Standard deviation of the Maxwellian distribution to sample
 * velocities from.
 * @param species_id Integer to assign to the "INTERNAL_STATE" property of each
 * particle.
 */

template <size_t ndim>
inline ParticleSet
uniform_cellwise_maxwellian(SYCLTargetSharedPtr sycl_target,
                            std::shared_ptr<PetscInterface::DMPlexInterface> mesh,
                            const ParticleSpec& particle_spec, const int npart_per_cell,
                            const REAL& weight, const REAL& std_dev,
                            const INT& species_id) {

  const int rank = sycl_target->comm_pair.rank_parent;

  std::mt19937 rng_pos(static_cast<std::mt19937::result_type>(52234234 + rank));
  std::mt19937 rng_vel(static_cast<std::mt19937::result_type>(52234231 + rank));

  std::vector<std::vector<double>> positions;
  std::vector<int> cell_ids;
  PetscInterface::uniform_within_dmplex_cells(mesh, npart_per_cell, positions, cell_ids,
                                              &rng_pos);

  const int N = static_cast<int>(cell_ids.size());

  auto velocities = NESO::Particles::normal_distribution(N, ndim, 0.0, std_dev, rng_vel);

  ParticleSet maxwellian(N, particle_spec);

  for (int px = 0; px < N; px++) {
    std::size_t pxu = static_cast<std::size_t>(px);
    for (int dimx = 0; dimx < static_cast<int>(ndim); dimx++) {
      std::size_t dimxu = static_cast<std::size_t>(dimx);
      maxwellian[Sym<REAL>("POSITION")][px][dimx] = positions.at(dimxu).at(pxu);
      maxwellian[Sym<REAL>("VELOCITY")][px][dimx] = velocities.at(dimxu).at(pxu);
    }
    maxwellian[Sym<INT>("CELL_ID")][px][0] = cell_ids.at(pxu);
    maxwellian[Sym<INT>("ID")][px][0] = px;
    maxwellian[Sym<REAL>("WEIGHT")][px][0] = weight;
    maxwellian[Sym<INT>("INTERNAL_STATE")][px][0] = species_id;
  }

  return maxwellian;
}

/**
 * @brief create an RNG kernel to use for sampling velocity distributions
 * in recombination and charge exchange.
 */

inline auto get_uniform_rng_kernel(SYCLTargetSharedPtr sycl_target, std::size_t n_samples,
                                   std::uint64_t root_seed = 141351) {

  const int rank = sycl_target->comm_pair.rank_parent;

  std::uint64_t seed = NESO::RNGToolkit::create_seeds(
      static_cast<std::size_t>(sycl_target->comm_pair.size_parent),
      static_cast<std::size_t>(rank), root_seed);

  auto rng_normal = NESO::RNGToolkit::create_rng<REAL>(
      NESO::RNGToolkit::Distribution::Uniform<REAL>{
          NESO::RNGToolkit::Distribution::next_value(0.0), 1.0},
      seed, sycl_target->device, static_cast<std::size_t>(sycl_target->device_index));

  // Create an interface between NESO-RNG-Toolkit and NESO-Particles KernelRNG
  auto rng_interface =
      make_rng_generation_function<GenericDeviceRNGGenerationFunction, REAL>(
          [=](REAL* d_ptr, const std::size_t num_samples) -> int {
            return rng_normal->get_samples(d_ptr, num_samples);
          });

  auto rng_kernel =
      host_atomic_block_kernel_rng<REAL>(rng_interface, static_cast<int>(n_samples));

  return rng_kernel;
}

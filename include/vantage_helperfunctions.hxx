#pragma once
#include "bout/bout.hxx"
#include <neso_particles.hpp>

using namespace NESO::Particles;

// get the number of cells owned by the local
// BOUT++ mesh, excluding guard cells
size_t get_num_cells_owned_bout_mesh(Mesh*& bout_mesh);

// get the cell volumes from the NESO-Particles mesh
// onto the same degrees of freedom owned by the local BOUT++ mesh
std::vector<REAL> get_cell_volumes_on_plasma_grid(
    DM& dm, std::vector<PetscInt>& kinetic_mesh_map,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh, Mesh*& bout_mesh);

// functions used for checks of the DMPlex

// Check that the x, y cell volumes of BOUT++ match the
// inferred quad cell volumes from the NESO-Particles mesh
void check_cell_volumes(std::vector<REAL> neso_cell_volumes_bmsh, Mesh*& bout_mesh,
                        Options& alloptions);

// Check that the Rxy, Zxy cell centres of BOUT++ match the
// inferred quad cell centres from the NESO-Particles mesh
void check_cell_centres(Options& alloptions, DM& dm,
                        std::vector<PetscInt>& kinetic_mesh_map,
                        std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        Mesh*& bout_mesh, BoutReal absolute_tolerance,
                        BoutReal relative_tolerance);

// Function to check mass conservation at the end of particle pushing
void check_mass_conservation(REAL total_mass_final, REAL total_mass_initial);

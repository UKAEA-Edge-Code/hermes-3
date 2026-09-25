#include "../include/vantage_helperfunctions.hxx"
#include "../include/component.hxx"
#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include <neso_particles.hpp>

using namespace NESO::Particles;

size_t get_num_cells_owned_bout_mesh(Mesh*& bout_mesh) {
  // local number of BOUT++ x cells, excluding guards
  const int Nx = bout_mesh->xend - bout_mesh->xstart + 1;
  // local number of BOUT++ y cells, excluding guards
  const int Ny = bout_mesh->yend - bout_mesh->ystart + 1;
  // Get the number of cells in the bout (plasma) mesh owned on this process, excluding guard cells
  const size_t num_cells_owned_bout_mesh = static_cast<size_t>(Nx * Ny);
  return num_cells_owned_bout_mesh;
}

std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0>
get_mesh_coupler_constant_weights(DM& dm, std::vector<PetscInt>& kinetic_mesh_map,
                                  Mesh*& bout_mesh, REAL backward_weight_0,
                                  REAL backward_weight_1) {
  const size_t num_cells_owned_bout_mesh = get_num_cells_owned_bout_mesh(bout_mesh);
  std::vector<std::vector<PetscInterface::DMPlexMeshCouplerDG0MapEntry>> coupler_map_0(
      static_cast<size_t>(num_cells_owned_bout_mesh));
  Field2D map_RZ_to_itriangle_0;
  Field2D map_RZ_to_itriangle_1;
  bout_mesh->get(map_RZ_to_itriangle_0, "map_RZ_to_itriangle_0");
  bout_mesh->get(map_RZ_to_itriangle_1, "map_RZ_to_itriangle_1");
  int icell = 0;
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      // lower triangle
      coupler_map_0.at(static_cast<size_t>(icell))
          .push_back(
              {kinetic_mesh_map.at(static_cast<size_t>(map_RZ_to_itriangle_0(ix, iy))),
               1.0, backward_weight_0});
      // upper triangle
      coupler_map_0.at(static_cast<size_t>(icell))
          .push_back(
              {kinetic_mesh_map.at(static_cast<size_t>(map_RZ_to_itriangle_1(ix, iy))),
               1.0, backward_weight_1});
      icell += 1;
    }
  }
  // object for transferring data between kinetic and bout mesh degree-of-freedom vectors
  std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0> mesh_coupler =
      std::make_shared<PetscInterface::DMPlexMeshCouplerDG0>(dm, coupler_map_0);
  return mesh_coupler;
}

std::vector<REAL> get_cell_volumes_on_plasma_grid(
    DM& dm, std::vector<PetscInt>& kinetic_mesh_map,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh, Mesh*& bout_mesh) {
  const size_t num_cells_owned_bout_mesh = get_num_cells_owned_bout_mesh(bout_mesh);
  // Get the number of cells in the kinetic (neutral) mesh owned on this process
  const size_t num_cells_owned_kinetic_mesh =
      static_cast<size_t>(neso_mesh->get_cell_count());
  // neso_mesh cell volumes on BOUT++ mesh indices
  std::vector<REAL> neso_cell_volumes_bmsh(num_cells_owned_bout_mesh);
  // the checks
  if (num_cells_owned_kinetic_mesh == num_cells_owned_bout_mesh) {
    // zero the compound index
    size_t ixy = 0;
    for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
      for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
        neso_cell_volumes_bmsh.at(ixy) =
            neso_mesh->dmh->get_cell_volume(static_cast<int>(ixy));
        ixy++;
      }
    }
  } else if (num_cells_owned_kinetic_mesh > num_cells_owned_bout_mesh) {
    // assume that this corresponds to the case where the BOUT++ mesh is decomposed
    // to triangles and there are also cells representing the region beyond the simulated plasma
    // -------------------------------------------
    // first, make a mesh_coupler_dg0 object with unit weights
    const std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0> mesh_coupler_unit_weight =
        get_mesh_coupler_constant_weights(dm, kinetic_mesh_map, bout_mesh, 1.0, 1.0);
    // obtain a list of kinetic mesh cell volumes
    std::vector<double> neso_cell_volumes_kmsh(num_cells_owned_kinetic_mesh);
    for (size_t ic = 0; ic < num_cells_owned_kinetic_mesh; ic++) {
      neso_cell_volumes_kmsh.at(ic) =
          neso_mesh->dmh->get_cell_volume(static_cast<int>(ic));
    }
    // move these cell volumes to the bout mesh
    mesh_coupler_unit_weight->backward_transfer(neso_cell_volumes_kmsh, 1,
                                                neso_cell_volumes_bmsh);
  }
  return neso_cell_volumes_bmsh;
}

std::vector<REAL> get_cell_vertices_on_plasma_grid(
    DM& dm, std::vector<PetscInt>& kinetic_mesh_map,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh, Mesh*& bout_mesh) {
  const size_t num_cells_owned_bout_mesh = get_num_cells_owned_bout_mesh(bout_mesh);
  // Get the number of cells in the kinetic (neutral) mesh owned on this process
  const size_t num_cells_owned_kinetic_mesh =
      static_cast<size_t>(neso_mesh->get_cell_count());
  // neso_mesh cell volumes on BOUT++ mesh indices
  const size_t nquad_vertices = 4;
  const size_t ntri_vertices = 3;
  const size_t ndim = 2; // number of position coordinates expected
  std::vector<REAL> quad_cell_vertices_bmsh(nquad_vertices * ndim
                                            * num_cells_owned_bout_mesh);
  // get the cell vertices in flattened vectors,
  // without attempting to respect anti-clockwise vertex ordering
  if (num_cells_owned_kinetic_mesh == num_cells_owned_bout_mesh) {
    std::vector<std::vector<REAL>> cell_vertices;
    // zero the compound index
    size_t ixy = 0;
    for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
      for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
        neso_mesh->dmh->get_cell_vertices(static_cast<PetscInt>(ixy), cell_vertices);
        for (size_t iv = 0; iv < nquad_vertices; iv++) {
          for (size_t idim = 0; idim < ndim; idim++) {
            const size_t jc = (ndim * ((nquad_vertices * ixy) + iv)) + idim;
            quad_cell_vertices_bmsh.at(jc) = cell_vertices.at(iv).at(idim);
          }
        }
        ixy++;
      }
    }
  } else if (num_cells_owned_kinetic_mesh > num_cells_owned_bout_mesh) {
    // assume that this corresponds to the case where the BOUT++ mesh is decomposed
    // to triangles and there are also cells representing the region beyond the simulated plasma
    // -------------------------------------------
    // we need to get the triangular cell coordinates from each upper and lower triangle
    // on to the local BOUT++ grid, then resolve which coordinates are unique to form
    // the coordinates for the quadrilateral cell which the pair of triangles represent
    // -------------------------------------------
    // first, make a mesh_coupler_dg0 object with unit weights from the lower triangle, and zero weight
    // for the upper triangle
    const std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0> mesh_coupler_0 =
        get_mesh_coupler_constant_weights(dm, kinetic_mesh_map, bout_mesh, 1.0, 0.0);
    // second, make a mesh_coupler_dg0 object with unit weights from the upper triangle, and zero weight
    // for the lower triangle
    const std::shared_ptr<PetscInterface::DMPlexMeshCouplerDG0> mesh_coupler_1 =
        get_mesh_coupler_constant_weights(dm, kinetic_mesh_map, bout_mesh, 0.0, 1.0);
    // obtain the cell coordinates for lower and upper triangles on the kinetic mesh
    std::vector<std::vector<REAL>> cell_vertices;
    std::vector<REAL> tri_cell_vertices_kmsh(ntri_vertices * ndim
                                             * num_cells_owned_kinetic_mesh);
    for (size_t ic = 0; ic < num_cells_owned_kinetic_mesh; ic++) {
      neso_mesh->dmh->get_cell_vertices(static_cast<PetscInt>(ic), cell_vertices);
      // fill in results to flattened vector
      for (size_t iv = 0; iv < ntri_vertices; iv++) {
        for (size_t idim = 0; idim < ndim; idim++) {
          const size_t jc = (ndim * ((ntri_vertices * ic) + iv)) + idim;
          tri_cell_vertices_kmsh.at(jc) = cell_vertices.at(iv).at(idim);
        }
      }
    }
    // transfer these results to vectors for the lower and upper triangles
    std::vector<REAL> tri_cell_vertices_0_bmsh(ntri_vertices * ndim
                                               * num_cells_owned_bout_mesh);
    std::vector<REAL> tri_cell_vertices_1_bmsh(ntri_vertices * ndim
                                               * num_cells_owned_bout_mesh);
    mesh_coupler_0->backward_transfer(tri_cell_vertices_kmsh, ntri_vertices * ndim,
                                      tri_cell_vertices_0_bmsh);
    mesh_coupler_1->backward_transfer(tri_cell_vertices_kmsh, ntri_vertices * ndim,
                                      tri_cell_vertices_1_bmsh);
    // fill in data for quad cell vertices
    // no requirement for the cell centre check to list in anti-clockwise order
    size_t ixy = 0;
    for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
      for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
        // first three vertices from lower triangle are definitely unqiue vertices for the quad
        // (though perhaps in an incorrect order)
        for (size_t iv = 0; iv < ntri_vertices; iv++) {
          for (size_t idim = 0; idim < ndim; idim++) {
            const size_t jc_quad = (ndim * ((nquad_vertices * ixy) + iv)) + idim;
            const size_t jc_tri = (ndim * ((ntri_vertices * ixy) + iv)) + idim;
            quad_cell_vertices_bmsh.at(jc_quad) = tri_cell_vertices_0_bmsh.at(jc_tri);
          }
        }
        // the final unique coordinate must be determined by checking for uniqueness
        const size_t ivquad = 3;
        const REAL atol = 1.0e-12;
        std::vector<bool> unique(ntri_vertices);
        for (size_t ivp = 0; ivp < ntri_vertices; ivp++) {
          const size_t jcp_tri = (ndim * ((ntri_vertices * ixy) + ivp));
          // initially presume that this index is unique
          unique.at(ivp) = true;
          for (size_t iv = 0; iv < ntri_vertices; iv++) {
            const size_t jc_tri = (ndim * ((ntri_vertices * ixy) + iv));
            REAL sumsqr = 0.0;
            // sum the squared lengths measuring the distance of this vertex from another
            for (size_t idim = 0; idim < ndim; idim++) {
              sumsqr += std::pow(tri_cell_vertices_0_bmsh.at(jc_tri + idim)
                                     - tri_cell_vertices_1_bmsh.at(jcp_tri + idim),
                                 2);
            }
            const REAL l2norm = std::sqrt(sumsqr);
            if (l2norm < atol) {
              unique.at(ivp) = false;
            }
          }
          if (unique.at(ivp)) {
            // this vertex has proved to be unique by not matching any other vertex
            for (size_t idim = 0; idim < ndim; idim++) {
              const size_t jc_quad = (ndim * ((nquad_vertices * ixy) + ivquad)) + idim;
              quad_cell_vertices_bmsh.at(jc_quad) =
                  tri_cell_vertices_1_bmsh.at(jcp_tri + idim);
            }
            // only one vertex can be unique
            break;
          }
        }
        ixy++;
      }
    }
  }
  return quad_cell_vertices_bmsh;
}

void check_cell_volumes(std::vector<REAL> neso_cell_volumes_bmsh, Mesh*& bout_mesh,
                        Options& alloptions) {
  Coordinates* coord = bout_mesh->getCoordinates();
  size_t ixy = 0;
  const REAL tolerance = 1.0e-12;
  const size_t num_cells_owned_bout_mesh = get_num_cells_owned_bout_mesh(bout_mesh);
  // Get the number of cells in the kinetic (neutral) mesh owned on this process
  ASSERT1(neso_cell_volumes_bmsh.size() == num_cells_owned_bout_mesh);
  // dimensional units
  const BoutReal meters = get<BoutReal>(alloptions["units"]["meters"]);
  const BoutReal meters_squared = meters * meters;
  const BoutReal meters_cubed = meters * meters * meters;
  // the checks of cell volumes
  // zero the compound index
  ixy = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      // Convert to SI: dx is m^2 T, J is m/T, dy is unitless, skip dz
      // so J * dx * dy = m^3, technically per radian toroidal angle due to missing dz
      const BoutReal bout_cell_area =
          coord->J(ix, iy) * coord->dx(ix, iy) * coord->dy(ix, iy) * meters_cubed;
      // neso_mesh is a 2D grid, needs m^2
      const REAL neso_cell_area = neso_cell_volumes_bmsh.at(ixy) * meters_squared;
      const bool volumes_match = (abs(bout_cell_area - neso_cell_area) < tolerance);
      // exit if we fail to find a match
      NESOASSERT(volumes_match,
                 fmt::format("BOUT++ mesh volume {} does not match NESO-Particles mesh "
                             "volume {} for ix = {} iy = {} \n Ignore this message by "
                             "setting [dmplex] test_dmplex_cell_volumes = false",
                             bout_cell_area, neso_cell_area, ix, iy));
      ixy++;
    }
  }
}

REAL cell_length(std::vector<std::vector<REAL>>& cell_vertices, std::size_t iv1,
                 std::size_t iv2, std::size_t iv3, std::size_t iv4) {
  const REAL Rlength2 =
      std::pow(0.5
                   * (cell_vertices.at(iv1).at(0) + cell_vertices.at(iv2).at(0)
                      - cell_vertices.at(iv3).at(0) - cell_vertices.at(iv4).at(0)),
               2.0);
  const REAL Zlength2 =
      std::pow(0.5
                   * (cell_vertices.at(iv1).at(1) + cell_vertices.at(iv2).at(1)
                      - cell_vertices.at(iv3).at(1) - cell_vertices.at(iv4).at(1)),
               2.0);
  REAL length = std::pow(Zlength2 + Rlength2, 0.5);
  return length;
}

void check_cell_centres(Options& alloptions, DM& dm,
                        std::vector<PetscInt>& kinetic_mesh_map,
                        std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        Mesh*& bout_mesh, BoutReal absolute_tolerance,
                        BoutReal relative_tolerance) {
  // get (R,Z) of cell centres in Hypnotoad grid
  Field2D Rxy;
  Field2D Zxy;
  bout_mesh->get(Rxy, "Rxy");
  bout_mesh->get(Zxy, "Zxy");

  BoutReal meters = get<BoutReal>(alloptions["units"]["meters"]);
  std::vector<REAL> neso_cell_vertices_plasma_grid =
      get_cell_vertices_on_plasma_grid(dm, kinetic_mesh_map, neso_mesh, bout_mesh);
  // number of vertices per quad
  const size_t nquad_vertices = 4;
  // expected dimensionality
  const size_t ndim = 2;
  // compare to cell centres calculated from cell corners
  std::vector<std::vector<REAL>> cell_vertices(nquad_vertices, std::vector<REAL>(ndim));

  PetscInt ixy = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      const REAL bout_Rxy = Rxy(ix, iy);
      const REAL bout_Zxy = Zxy(ix, iy);

      // fill in the vertices from the flattened vector
      for (std::size_t iv = 0; iv < nquad_vertices; iv++) {
        for (size_t idim = 0; idim < ndim; idim++) {
          const size_t jc_quad =
              (ndim * ((nquad_vertices * static_cast<size_t>(ixy)) + iv)) + idim;
          cell_vertices.at(iv).at(idim) = neso_cell_vertices_plasma_grid.at(jc_quad);
        }
      }
      REAL neso_Rxy = 0.0;
      REAL neso_Zxy = 0.0;
      for (std::size_t iv = 0; iv < nquad_vertices; iv++) {
        // DMPlex is stored in normalised units, need conversion to [m]
        neso_Rxy += cell_vertices.at(iv).at(0) * meters;
        neso_Zxy += cell_vertices.at(iv).at(1) * meters;
      }
      neso_Rxy /= 4.0;
      neso_Zxy /= 4.0;
      // get lengths of cell across the two dimensions
      const REAL cell_length_a = cell_length(cell_vertices, 0, 1, 2, 3) * meters;
      const REAL cell_length_b = cell_length(cell_vertices, 0, 3, 2, 1) * meters;
      const REAL min_cell_length = std::min(cell_length_a, cell_length_b);
      // we compare the difference in cell centres to the absolute tolerance and
      // the relative tolerance formed by comparing to the smallest length across the cell
      const REAL tolerance = absolute_tolerance + min_cell_length * relative_tolerance;
      const bool centres_match = (abs(neso_Rxy - bout_Rxy) < tolerance)
                                 && (abs(neso_Zxy - bout_Zxy) < tolerance);
      // exit if we fail to find a match
      NESOASSERT(
          centres_match,
          fmt::format("Hypnotoad/BOUT++ cell centre (R, Z) ({}, {}) does not match "
                      "NESO-Particles mesh inferred quad "
                      "cell centre ({}, {}) for ix = {} iy = {} \n"
                      "The cell height and width are {} {} \n"
                      "The displacements in R and Z are {} {} \n"
                      "Ignore this message by "
                      "setting [dmplex] test_dmplex_cell_centres = false\n Relax the "
                      "tolerance used in this check by increasing\n"
                      "[dmplex] dmplex_cell_centre_absolute_tolerance = {}\n"
                      "[dmplex] dmplex_cell_centre_relative_tolerance = {}",
                      bout_Rxy, bout_Zxy, neso_Rxy, neso_Zxy, ix, iy, cell_length_a,
                      cell_length_b, abs(neso_Rxy - bout_Rxy), abs(neso_Zxy - bout_Zxy),
                      absolute_tolerance, relative_tolerance));
      ixy++;
    }
  }
}

void check_mass_conservation(REAL total_mass_final, REAL total_mass_initial) {
  REAL rtol = 1.0e-13;
  REAL mass_conserved =
      (abs(total_mass_final - total_mass_initial) < rtol * total_mass_initial);
  // exit if we fail to find conservation
  NESOASSERT(mass_conserved,
             fmt::format("Initial total mass {} does not match "
                         "final total mass {} \n Ignore this message by "
                         "setting [vantage] test_mass_conservation = false",
                         total_mass_initial, total_mass_final));
}

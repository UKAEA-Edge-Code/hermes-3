#include "../include/vantage_dmplex.hxx"
#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include <bout/assert.hxx>
#include <bout/field_factory.hxx>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fmt/core.h>
#include <fmt/format.h>
#include <neso_particles.hpp>
#include <neso_particles/compute_target.hpp>
#include <neso_particles/containers/cell_data.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include <neso_particles/typedefs.hpp>
#include <netcdf>
#include <petscsystypes.h>
#include <petscviewerhdf5.h>
#include <vector>

#ifndef NESO_PARTICLES_PETSC
static_assert(false, "NESO-Particles was installed without PETSc support.");
#else

using namespace NESO::Particles;

template <typename T, typename U>
inline void ASSERT_EQ(T t, U u) {
  NESOASSERT(t == u, "A check failed.");
}

void collect_unique_points(std::vector<double>& global_Z_vertices_buffer,
                           std::vector<double>& global_R_vertices_buffer,
                           size_t& N_unique, const double& tolerance,
                           std::vector<double>& global_Z_hypnotoad_vertices,
                           std::vector<double>& global_R_hypnotoad_vertices) {
  bool unique;
  size_t N_nonunique_vertices = global_Z_hypnotoad_vertices.size();
  for (size_t iv = 0; iv < N_nonunique_vertices; iv++) {
    // assume the point
    // iv = (global_R_hypnotoad_vertices.at(iv),global_Z_hypnotoad_vertices.at(iv))
    // is unique
    unique = true;
    // check if the point is unique, by comparing the the existing N_unique points
    for (size_t iunique = 0; iunique < N_unique; iunique++) {
      if (std::abs(global_Z_hypnotoad_vertices.at(iv)
                   - global_Z_vertices_buffer.at(iunique))
              < tolerance
          && std::abs(global_R_hypnotoad_vertices.at(iv)
                      - global_R_vertices_buffer.at(iunique))
                 < tolerance) {
        unique = false;
        // we have determined that the point is not unique
      }
    }
    if (unique) {
      // add the point and increment N_unique
      global_Z_vertices_buffer.at(N_unique) = global_Z_hypnotoad_vertices.at(iv);
      global_R_vertices_buffer.at(N_unique) = global_R_hypnotoad_vertices.at(iv);
      N_unique++;
    }
  }
}

void RZ_to_ivertex_vector(Field2D& ivertex_corners,
                          std::vector<double>& global_Z_vertices,
                          std::vector<double>& global_R_vertices, const double& tolerance,
                          Mesh*& bout_mesh, Field2D& Rxy_corners, Field2D& Zxy_corners) {
  size_t Nvertex = global_Z_vertices.size();
  bool index_found;
  // loop over points that are not guard cells
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      index_found = false;
      for (size_t iv = 0; iv < Nvertex; iv++) {
        if (std::abs(global_Z_vertices.at(iv) - Zxy_corners(ix, iy)) < tolerance
            && std::abs(global_R_vertices.at(iv) - Rxy_corners(ix, iy)) < tolerance) {
          ivertex_corners(ix, iy) = static_cast<BoutReal>(iv);
          index_found = true;
          // we have matched an iv global index to a (R,Z) from Hypnotoad
        }
      }
      // exit if we fail to find a match
      NESOASSERT(index_found,
                 fmt::format("ivertex not found for ix = {} iy = {}", ix, iy));
    }
  }
}

VantageBasicMeshData cells_definition_from_RZ_ivertex(
    Mesh*& bout_mesh, Field2D& Rxy_lower_left_corners, Field2D& Rxy_lower_right_corners,
    Field2D& Rxy_upper_right_corners, Field2D& Rxy_upper_left_corners,
    Field2D& Zxy_lower_left_corners, Field2D& Zxy_lower_right_corners,
    Field2D& Zxy_upper_right_corners, Field2D& Zxy_upper_left_corners,
    std::vector<double>& global_R_vertices, std::vector<double>& global_Z_vertices,
    const BoutReal dmplex_vertex_tolerance) {
  ASSERT1(global_R_vertices.size() == global_Z_vertices.size());
  // ivertex arrays made in cxx, initialise with -1 index
  Field2D ivertex_lower_left_corners{-1, bout_mesh};
  Field2D ivertex_lower_right_corners{-1, bout_mesh};
  Field2D ivertex_upper_right_corners{-1, bout_mesh};
  Field2D ivertex_upper_left_corners{-1, bout_mesh};
  // now fill ivertex_corners arrays
  // these arrays identify a given (R,Z) location with one of the vertices in the global list
  RZ_to_ivertex_vector(ivertex_lower_left_corners, global_Z_vertices, global_R_vertices,
                       dmplex_vertex_tolerance, bout_mesh, Rxy_lower_left_corners,
                       Zxy_lower_left_corners);
  RZ_to_ivertex_vector(ivertex_lower_right_corners, global_Z_vertices, global_R_vertices,
                       dmplex_vertex_tolerance, bout_mesh, Rxy_lower_right_corners,
                       Zxy_lower_right_corners);
  RZ_to_ivertex_vector(ivertex_upper_right_corners, global_Z_vertices, global_R_vertices,
                       dmplex_vertex_tolerance, bout_mesh, Rxy_upper_right_corners,
                       Zxy_upper_right_corners);
  RZ_to_ivertex_vector(ivertex_upper_left_corners, global_Z_vertices, global_R_vertices,
                       dmplex_vertex_tolerance, bout_mesh, Rxy_upper_left_corners,
                       Zxy_upper_left_corners);
  // use the gloabl list of vertices, and their identification in the Field2D ivertex arrays
  // to construct anticlockwise listed quads, and split these into anticlockwise listed
  // triangular cells
  std::vector<PetscReal> Z_vertices(4);
  std::vector<PetscReal> R_vertices(4);
  std::vector<PetscReal> theta_vertices(4);
  std::vector<PetscInt> i_vertices(4);
  std::vector<size_t> sort_indices(4);
  PetscReal ZZ;
  PetscReal ZZmid;
  PetscReal RR;
  PetscReal RRmid;
  // local number of x cells, excluding guards
  const int Nx = bout_mesh->xend - bout_mesh->xstart + 1;
  // local number of y cells, excluding guards
  const int Ny = bout_mesh->yend - bout_mesh->ystart + 1;
  const PetscInt num_quad_cells_owned = Nx * Ny;
  const int nranks = BoutComm::size();
  const int irank = BoutComm::rank();
  // number of quad cells in Bout mesh is number of ranks times the number of cells per rank
  // number of triangular cells inferred from this is 2 * nranks * num_cells_owned
  // number of vertices per cell is 3
  const int ntri_vertices = 3;
  const int ntriangles_per_rank = 2 * num_quad_cells_owned;
  const int nvertices_per_rank = ntri_vertices * ntriangles_per_rank;
  const size_t shift = static_cast<size_t>(irank * nvertices_per_rank);
  // the two triangles [0, 1, 2], [0, 2, 3] -> anticlockwise if [0, 1, 2, 3] anticlockwise
  const std::vector<size_t> itriangle_0{0, 1, 2};
  const std::vector<size_t> itriangle_1{0, 2, 3};
  std::vector<PetscInt> cells_local(static_cast<size_t>(nvertices_per_rank * nranks),
                                    0.0);
  // maps from (R,Z) to global triangle index
  Field2D map_RZ_to_itriangle_0{-1, bout_mesh};
  Field2D map_RZ_to_itriangle_1{-1, bout_mesh};
  // We are careful to list the vertices in counter clock-wise order.
  int ixy = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      // collect data from Hypnotoad arrays
      i_vertices[0] = static_cast<int>(std::lround(ivertex_lower_left_corners(ix, iy)));
      i_vertices[1] = static_cast<int>(std::lround(ivertex_lower_right_corners(ix, iy)));
      i_vertices[2] = static_cast<int>(std::lround(ivertex_upper_right_corners(ix, iy)));
      i_vertices[3] = static_cast<int>(std::lround(ivertex_upper_left_corners(ix, iy)));

      R_vertices[0] = Rxy_lower_left_corners(ix, iy);
      R_vertices[1] = Rxy_lower_right_corners(ix, iy);
      R_vertices[2] = Rxy_upper_right_corners(ix, iy);
      R_vertices[3] = Rxy_upper_left_corners(ix, iy);

      Z_vertices[0] = Zxy_lower_left_corners(ix, iy);
      Z_vertices[1] = Zxy_lower_right_corners(ix, iy);
      Z_vertices[2] = Zxy_upper_right_corners(ix, iy);
      Z_vertices[3] = Zxy_upper_left_corners(ix, iy);
      // get the mean values of R, Z
      RRmid = 0.0;
      ZZmid = 0.0;
      for (size_t iv = 0; iv < 4; iv++) {
        RRmid += R_vertices[iv];
        ZZmid += Z_vertices[iv];
      }
      RRmid /= 4.0;
      ZZmid /= 4.0;
      // get the angle subtended from the centre of the cell to each vertex
      for (size_t iv = 0; iv < 4; iv++) {
        RR = R_vertices[iv] - RRmid;
        ZZ = Z_vertices[iv] - ZZmid;
        theta_vertices[iv] = std::atan2(ZZ, RR);
      }
      // get the indices that sort the vertices in ascending order of theta
      for (size_t iv = 0; iv < 4; ++iv) {
        sort_indices[iv] = iv;
      }
      std::sort(sort_indices.begin(), sort_indices.end(),
                [&theta_vertices](size_t i, size_t j) {
                  return theta_vertices[i] < theta_vertices[j];
                });

      // fill cells using the sorted indices
      // noting that quad cell vertices defined anticlockwise with indices [0, 1, 2 ,3]
      // mean that the two triangular cell vertices are defined anticlockwise as [0, 1, 2], [0, 2, 3]
      const size_t shift_inner = shift + (static_cast<size_t>(ntri_vertices * 2 * ixy));
      for (size_t iv = 0; iv < static_cast<size_t>(ntri_vertices); ++iv) {
        const size_t itri_0 = itriangle_0.at(iv);
        const size_t itri_1 = itriangle_1.at(iv);
        // assign cell definition for cells on the local process
        cells_local.at(shift_inner + iv) = (i_vertices[sort_indices[itri_0]]);
        cells_local.at(shift_inner + static_cast<size_t>(ntri_vertices) + iv) =
            (i_vertices[sort_indices[itri_1]]);
      }
      map_RZ_to_itriangle_0(ix, iy) = 2 * ixy + ntriangles_per_rank * irank;
      map_RZ_to_itriangle_1(ix, iy) = 2 * ixy + 1 + ntriangles_per_rank * irank;
      ixy++;
    }
  }
  // use MPIAllreduce to get the global cell definitions on all ranks
  std::vector<PetscInt> cells(cells_local.size(), 0.0);
  MPICHK(MPI_Allreduce(cells_local.data(), cells.data(), static_cast<int>(cells.size()),
                       MPI_INT, MPI_SUM, BoutComm::get()));
  // make a flattened vector of the global vertices lists
  std::vector<double> vertices(2 * global_R_vertices.size());
  for (size_t iv = 0; iv < global_R_vertices.size(); iv++) {
    vertices.at((2 * iv) + 0) = global_R_vertices.at(iv);
    vertices.at((2 * iv) + 1) = global_Z_vertices.at(iv);
  }
  return VantageBasicMeshData{
      map_RZ_to_itriangle_0, map_RZ_to_itriangle_1,
      VerticesData{vertices, global_R_vertices.size(), 2},
      TrianglesDefinitionData{cells, static_cast<size_t>(nranks * ntriangles_per_rank),
                              ntri_vertices}};
}

void write_dmplex_to_file(DM& dm, const std::string& dmplex_name,
                          const std::string& dmplex_h5_filename) {
  // save a HDF5 file containing the DM for diagnostics
  PetscViewer viewer;
  // Set a name for the DMPlex object (important for HDF5)
  PetscObjectSetName(reinterpret_cast<PetscObject>(dm), dmplex_name.c_str());
  // Create an HDF5 viewer
  PetscViewerHDF5Open(BoutComm::get(), dmplex_h5_filename.c_str(), FILE_MODE_WRITE,
                      &viewer);
  // Set viewer format to PETSC_VIEWER_HDF5_PETSC for compatibility
  PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_PETSC);
  // Save the DMPlex to the HDF5 file
  DMView(dm, viewer);
  // Clean up
  PetscViewerDestroy(&viewer);
  output << "Finished DMPlex diagnostic \n";
}

VantageBasicMeshData kinetic_mesh_data_from_netcdf(Mesh* bout_mesh) {
  Field2D map_RZ_to_itriangle_0;
  Field2D map_RZ_to_itriangle_1;
  const int read_status_itri0 =
      bout_mesh->get(map_RZ_to_itriangle_0, "map_RZ_to_itriangle_0");
  ASSERT1(read_status_itri0 == 0) // check map read successfully from file
  const int read_status_itri1 =
      bout_mesh->get(map_RZ_to_itriangle_1, "map_RZ_to_itriangle_1");
  ASSERT1(read_status_itri1 == 0) // check map read successfully from file
  // get data that defines triangular cells
  const VerticesData vertices_data = get_triangle_vertices();
  const TrianglesDefinitionData cell_definition = get_triangle_cell_definition();
  return VantageBasicMeshData{map_RZ_to_itriangle_0, map_RZ_to_itriangle_1, vertices_data,
                              cell_definition};
}

VantageBasicMeshData kinetic_mesh_data_from_Bout_mesh(Mesh* bout_mesh,
                                                      Options& mesh_options) {

  // DMPlex vertex distance tolerance for duplicate Hypnotoad vertices
  const BoutReal dmplex_vertex_tolerance =
      mesh_options["dmplex_vertex_tolerance"]
          .doc("Tolerance for determining duplicate vertices when creating DMPlex from "
               "BOUT++ mesh.")
          .withDefault(1.0e-8);

  Field2D Rxy_lower_left_corners;
  Field2D Rxy_lower_right_corners;
  Field2D Rxy_upper_right_corners;
  Field2D Rxy_upper_left_corners;
  Field2D Zxy_lower_left_corners;
  Field2D Zxy_lower_right_corners;
  Field2D Zxy_upper_right_corners;
  Field2D Zxy_upper_left_corners;
  const int read_status_Rxy = bout_mesh->get(Rxy_lower_left_corners, "Rxy_corners");
  ASSERT1(read_status_Rxy == 0);
  const int read_status_Rxy_lr =
      bout_mesh->get(Rxy_lower_right_corners, "Rxy_lower_right_corners");
  ASSERT1(read_status_Rxy_lr == 0);
  const int read_status_Rxy_ur =
      bout_mesh->get(Rxy_upper_right_corners, "Rxy_upper_right_corners");
  ASSERT1(read_status_Rxy_ur == 0);
  const int read_status_Rxy_ul =
      bout_mesh->get(Rxy_upper_left_corners, "Rxy_upper_left_corners");
  ASSERT1(read_status_Rxy_ul == 0);
  const int read_status_Zxy = bout_mesh->get(Zxy_lower_left_corners, "Zxy_corners");
  ASSERT1(read_status_Zxy == 0);
  const int read_status_Zxy_lr =
      bout_mesh->get(Zxy_lower_right_corners, "Zxy_lower_right_corners");
  ASSERT1(read_status_Zxy_lr == 0);
  const int read_status_Zxy_ur =
      bout_mesh->get(Zxy_upper_right_corners, "Zxy_upper_right_corners");
  ASSERT1(read_status_Zxy_ur == 0);
  const int read_status_Zxy_ul =
      bout_mesh->get(Zxy_upper_left_corners, "Zxy_upper_left_corners");
  ASSERT1(read_status_Zxy_ul == 0);
  Field2D ivertex_lower_left_corners;
  Field2D ivertex_lower_right_corners;
  Field2D ivertex_upper_right_corners;
  Field2D ivertex_upper_left_corners;
  // local number of x cells, excluding guards
  const int Nx = bout_mesh->xend - bout_mesh->xstart + 1;
  // local number of y cells, excluding guards
  const int Ny = bout_mesh->yend - bout_mesh->ystart + 1;

  const int mpi_size = BoutComm::size();
  const int mpi_rank = BoutComm::rank();
  // global number of physical nonunique vertices stored in hypnotoad datasets
  const size_t N_nonunique_vertices = static_cast<size_t>(mpi_size * Nx * Ny);
  // arrays to fill with local data
  std::vector<double> local_Z_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_Z_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_Z_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_Z_upper_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_upper_left_vertices(N_nonunique_vertices, 0.0);
  // arrays to receive the summed data across ranks
  std::vector<double> global_Z_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_Z_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_Z_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_Z_upper_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_upper_left_vertices(N_nonunique_vertices, 0.0);
  // fill these vectors with vertex values from the local rank
  // at indices determined by the local rank
  size_t icxy = static_cast<size_t>(Nx * Ny * mpi_rank);
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      local_R_lower_left_vertices.at(icxy) = Rxy_lower_left_corners(ix, iy);
      local_Z_lower_left_vertices.at(icxy) = Zxy_lower_left_corners(ix, iy);
      local_R_lower_right_vertices.at(icxy) = Rxy_lower_right_corners(ix, iy);
      local_Z_lower_right_vertices.at(icxy) = Zxy_lower_right_corners(ix, iy);
      local_R_upper_right_vertices.at(icxy) = Rxy_upper_right_corners(ix, iy);
      local_Z_upper_right_vertices.at(icxy) = Zxy_upper_right_corners(ix, iy);
      local_R_upper_left_vertices.at(icxy) = Rxy_upper_left_corners(ix, iy);
      local_Z_upper_left_vertices.at(icxy) = Zxy_upper_left_corners(ix, iy);
      icxy++;
    }
  }
  // Perform Allreduce (sum) to get knowledge of vertices to all ranks
  MPICHK(MPI_Allreduce(
      local_R_lower_left_vertices.data(), global_R_lower_left_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  MPICHK(MPI_Allreduce(
      local_Z_lower_left_vertices.data(), global_Z_lower_left_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  MPICHK(MPI_Allreduce(
      local_R_lower_right_vertices.data(), global_R_lower_right_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  MPICHK(MPI_Allreduce(
      local_Z_lower_right_vertices.data(), global_Z_lower_right_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  MPICHK(MPI_Allreduce(
      local_R_upper_right_vertices.data(), global_R_upper_right_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  MPICHK(MPI_Allreduce(
      local_Z_upper_right_vertices.data(), global_Z_upper_right_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  MPICHK(MPI_Allreduce(
      local_R_upper_left_vertices.data(), global_R_upper_left_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  MPICHK(MPI_Allreduce(
      local_Z_upper_left_vertices.data(), global_Z_upper_left_vertices.data(),
      static_cast<int>(N_nonunique_vertices), MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  // if (mpi_rank == 0) {
  //     std::cout << "Result of Allreduce (sum): ";
  //     for (double val : global_R_lower_left_vertices) {
  //         std::cout << val << " ";
  //     }
  //     std::cout << std::endl;
  //     std::cout << "N_nonunique_vertices=" << N_nonunique_vertices << std::endl;
  // }
  // Now dynamically determine a list of unique vertex points
  // constant to give us a vector that can definitely contain all points in the global
  // lists
  const size_t N_global_nonunique_vertices = static_cast<size_t>(4 * mpi_size * Nx * Ny);
  std::vector<double> global_Z_vertices_buffer(N_global_nonunique_vertices, 0.0);
  std::vector<double> global_R_vertices_buffer(N_global_nonunique_vertices, 0.0);
  // fill the buffer vectors, checking each time if the point is unique
  // first point, outside loop
  global_Z_vertices_buffer.at(0) = global_Z_lower_left_vertices.at(0);
  global_R_vertices_buffer.at(0) = global_R_lower_left_vertices.at(0);
  size_t N_unique = 1; // we have one unique point in the buffer
  // loop over lower left vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        dmplex_vertex_tolerance, global_Z_lower_left_vertices,
                        global_R_lower_left_vertices);
  // loop over lower right vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        dmplex_vertex_tolerance, global_Z_lower_right_vertices,
                        global_R_lower_right_vertices);
  // loop over upper right vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        dmplex_vertex_tolerance, global_Z_upper_right_vertices,
                        global_R_upper_right_vertices);
  // loop over upper left vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        dmplex_vertex_tolerance, global_Z_upper_left_vertices,
                        global_R_upper_left_vertices);
  // now make a vector of the size N_unique and fill from the buffer
  std::vector<double> global_Z_vertices(N_unique, 0.0);
  std::vector<double> global_R_vertices(N_unique, 0.0);
  for (size_t iv = 0; iv < N_unique; iv++) {
    global_Z_vertices.at(iv) = global_Z_vertices_buffer.at(iv);
    global_R_vertices.at(iv) = global_R_vertices_buffer.at(iv);
  }
  if (mpi_rank == 0) {
    std::cout << "Result of vertex collection: ";
    // for (int iv=0; iv<N_unique; iv++) {
    //     std::cout << "(" << global_R_vertices.at(iv) << ", " <<
    //     global_Z_vertices.at(iv) << ") ";
    // }
    // std::cout << std::endl;
    std::cout << "N_unique=" << N_unique << std::endl;
  }

  return cells_definition_from_RZ_ivertex(
      bout_mesh, Rxy_lower_left_corners, Rxy_lower_right_corners, Rxy_upper_right_corners,
      Rxy_upper_left_corners, Zxy_lower_left_corners, Zxy_lower_right_corners,
      Zxy_upper_right_corners, Zxy_upper_left_corners, global_R_vertices,
      global_Z_vertices, dmplex_vertex_tolerance);
}

void create_dmplex_in_serial(VantageBasicMeshData& basic_mesh_data, DM& dm) {
  // First we setup the integers for the topology of the mesh (using data in serial only).
  const PetscInt num_cells_owned =
      static_cast<PetscInt>(basic_mesh_data.cell_definition.ntriangles_global);
  const PetscInt num_vertices_owned =
      static_cast<PetscInt>(basic_mesh_data.vertices_data.nvertices_global);
  const std::vector<PetscInt> cells = basic_mesh_data.cell_definition.tri_cell_vertices;
  const std::vector<PetscScalar> vertex_coords = basic_mesh_data.vertices_data.vertices;
  const PetscInt ndim = static_cast<PetscInt>(basic_mesh_data.vertices_data.ndim);
  const PetscInt ncorners =
      static_cast<PetscInt>(basic_mesh_data.cell_definition.ncorners);
  // create the DMPlex in serial
  PETSCCHK(DMPlexCreateFromCellListPetsc(BoutComm::get(), ndim, num_cells_owned,
                                         num_vertices_owned, ncorners, PETSC_TRUE,
                                         cells.data(), ndim, vertex_coords.data(), &dm));
}

VerticesData get_triangle_vertices() {
  // read data from netcdf for global vertices in mesh
  // Open the NetCDF file in read-only mode
  const std::string filename = Options::root()["mesh"]["file"];
  netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

  // Get the vertices variable
  // vertices is a global list of vertex coordinates
  // std::string varName = "vertices";
  netCDF::NcVar dataVar_vertices = dataFile.getVar("vertices");
  NESOASSERT(!dataVar_vertices.isNull(), "vertices not found in file.");
  std::vector<netCDF::NcDim> dims_vertices = dataVar_vertices.getDims();
  size_t nvertices = dims_vertices[0].getSize();
  size_t ncomp = dims_vertices[1].getSize();
  // Read the data into a vector
  std::vector<REAL> vertices(nvertices * ncomp);
  dataVar_vertices.getVar(vertices.data());

  // close the netcdf file
  dataFile.close();

  return VerticesData{vertices, nvertices, ncomp};
}

TrianglesDefinitionData get_triangle_cell_definition() {
  // read data from netcdf for global vertices in mesh
  // Open the NetCDF file in read-only mode
  const std::string filename = Options::root()["mesh"]["file"];
  netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

  // Get the tri_cell_vertices variable
  // a list of integers defining each triangular cell
  // in terms of indices that
  // index the "vertices" list loaded above
  // std::string varName_tri_cell = "tri_cell_vertices";
  netCDF::NcVar dataVar_tri_cell = dataFile.getVar("tri_cell_vertices");
  NESOASSERT(!dataVar_tri_cell.isNull(), "tri_cell_vertices not found in file.");
  std::vector<netCDF::NcDim> dims_tri_cell = dataVar_tri_cell.getDims();
  size_t ntriangle = dims_tri_cell[0].getSize();
  size_t ntricorners = dims_tri_cell[1].getSize();
  // Read the data into a vector
  std::vector<int> tri_cell_vertices(ntriangle * ntricorners);
  dataVar_tri_cell.getVar(tri_cell_vertices.data());
  // close the netcdf file
  dataFile.close();
  // std::cout << "tri_cell_verticies" << "\n";
  // for (size_t it=0; it < ntriangle; it++){
  //   std::cout << fmt::format("local_cell.at({}): ",it);
  //   for (size_t iv=0; iv < 3; iv++){
  //     std::cout << " " << tri_cell_vertices.at((it*3) + iv) << ", ";
  //   }
  //   std::cout << "\n ";
  // }
  return TrianglesDefinitionData{tri_cell_vertices, ntriangle, ntricorners};
}

REAL get_triangle_area(size_t itriangle, const std::vector<REAL>& vertices,
                       const std::vector<int>& tri_cell_vertices) {
  // compute the area for this triangle
  // use result of vector product for area
  // A = 1/2 | u x v |
  // where u and v are vectors defining two sides of the triangle

  // three vertices per triangle
  const size_t ntri = 3;
  std::vector<int> local_cell(ntri);
  // obtain the global vertex integers which define the local triangular cell
  for (size_t iv = 0; iv < local_cell.size(); iv++) {
    local_cell.at(iv) = tri_cell_vertices.at((ntri * itriangle) + iv);
    // std::cout << fmt::format("local_cell.at({}): ",iv) << local_cell.at(iv) << '\n';
  }
  // std::cout << "local_cell: " << local_cell.data() << '\n';
  // expect two vector components per vertex, mesh is 2D
  const size_t ncomp = 2;
  std::vector<double> local_vertices(ntri * ncomp);
  for (size_t iv = 0; iv < local_cell.size(); iv++) {
    for (size_t ic = 0; ic < ncomp; ic++) {
      const size_t jc = (iv * ncomp) + ic;
      local_vertices.at(jc) =
          vertices.at((static_cast<size_t>(local_cell.at(iv)) * ncomp) + ic);
      // std::cout << fmt::format("local_vertices.at({}): ",jc) << local_vertices.at(jc) << '\n';
    }
  }
  const size_t iv0 = 0;
  const size_t iv1 = 1;
  const size_t iv2 = 2;
  const REAL ux = local_vertices.at(iv1 * ncomp) - local_vertices.at(iv0);
  const REAL uy = local_vertices.at((iv1 * ncomp) + 1) - local_vertices.at(iv0 + 1);
  const REAL vx = local_vertices.at(iv2 * ncomp) - local_vertices.at(iv0);
  const REAL vy = local_vertices.at((iv2 * ncomp) + 1) - local_vertices.at(iv0 + 1);
  const REAL area = 0.5 * std::abs((ux * vy) - (uy * vx));
  // std::cout << "area: " << area << '\n';
  return area;
}

#endif

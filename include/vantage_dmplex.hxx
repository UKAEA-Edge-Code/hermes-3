#pragma once
#include <cstddef>
#include <string>
#ifndef VANTAGE_DMPLEX_H
#define VANTAGE_DMPLEX_H

#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include <bout/field_factory.hxx>
#include <cmath>
#include <fmt/core.h>
#include <neso_particles.hpp>
#include <neso_particles/compute_target.hpp>
#include <neso_particles/containers/cell_data.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include <neso_particles/typedefs.hpp>
#include <netcdf>
#include <petscsystypes.h>
#include <petscviewerhdf5.h>

#ifndef NESO_PARTICLES_PETSC
static_assert(false, "NESO-Particles was installed without PETSc support.");
#else

using namespace NESO::Particles;

// data recording the global vertex coordinates for the DMPlex mesh
struct VerticesData {
  // global, flattened list of vertices coordinates in the mesh
  std::vector<REAL> vertices;
  // total number of vertices in the mesh
  size_t nvertices_global;
  // number of dimensions in the mesh
  size_t ndim;
};

// data recording the definition of triangular DMPlex cells in the mesh
struct TrianglesDefinitionData {
  // flattened vector of integers defining the triangular cells in the mesh, with the indices defined by "vertices" above
  std::vector<int> tri_cell_vertices;
  // total number of triangular cells in the mesh
  size_t ntriangles_global;
  // number of corners in each cell in the mesh
  size_t ncorners;
};

// auxiliary information needed to create the
// mesh coupler object and perform other initialisation steps
// with the NESO-Particles kinetic mesh
struct VantageBasicMeshData {
  // map from R,Z (x,y) to the triangle "0" in a given BOUT++ cell
  Field2D map_RZ_to_itriangle_0;
  // map from R,Z (x,y) to the triangle "1" in a given BOUT++ cell
  Field2D map_RZ_to_itriangle_1;
  // vertex coordinates in the mesh
  VerticesData vertices_data;
  // data defining the triangular cells in the mesh
  TrianglesDefinitionData cell_definition;
};

VantageBasicMeshData cells_definition_from_RZ_ivertex(
    Mesh*& bout_mesh, Field2D& Rxy_lower_left_corners, Field2D& Rxy_lower_right_corners,
    Field2D& Rxy_upper_right_corners, Field2D& Rxy_upper_left_corners,
    Field2D& Zxy_lower_left_corners, Field2D& Zxy_lower_right_corners,
    Field2D& Zxy_upper_right_corners, Field2D& Zxy_upper_left_corners,
    std::vector<double>& global_R_vertices, std::vector<double>& global_Z_vertices,
    const BoutReal dmplex_vertex_tolerance);

VantageBasicMeshData create_dmplex_from_Bout_mesh(DM& dm, Mesh* bout_mesh,
                                                  Options& mesh_options);

VantageBasicMeshData create_dmplex_from_GMSH_msh(DM& dm, Mesh* bout_mesh,
                                                 std::string msh_file);

void write_dmplex_to_file(DM dm, std::string dmplex_name, std::string dmplex_h5_filename);

BoutReal get_triangle_area(size_t itriangle, const std::vector<double>& vertices,
                           const std::vector<int>& tri_cell_vertices);

VerticesData get_triangle_vertices();

TrianglesDefinitionData get_triangle_cell_definition();

#endif

#endif // VANTAGE_DMPLEX_H

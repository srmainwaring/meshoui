// ==========================================================================
// Helios
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "meshoui/meshoui.h"

using namespace meshoui;

// Test the two methods for loading a mesh.

int main() {

  // Two methods for loading: using vertices and faces (input_load = true) or using an input mesh file (input_load = false).

  // Building a box.
  std::vector<meshoui::Vector3d> vertices;
  vertices.emplace_back(-1, -1, 1);
  vertices.emplace_back(1, -1, 1);
  vertices.emplace_back(1, 1, 1);
  vertices.emplace_back(-1, 1, 1);
  vertices.emplace_back(-1, -1, -1);
  vertices.emplace_back(1, -1, -1);
  vertices.emplace_back(1, 1, -1);
  vertices.emplace_back(-1, 1, -1);
  std::vector<Eigen::VectorXi> faces;
  faces.emplace_back(Eigen::Vector3i(0, 1, 2));
  faces.emplace_back(Eigen::Vector3i(2, 3, 0));
  faces.emplace_back(Eigen::Vector3i(0, 4, 1));
  faces.emplace_back(Eigen::Vector3i(1, 4, 5));
  faces.emplace_back(Eigen::Vector3i(1, 5, 2));
  faces.emplace_back(Eigen::Vector3i(2, 5, 6));
  faces.emplace_back(Eigen::Vector3i(2, 6, 3));
  faces.emplace_back(Eigen::Vector3i(3, 6, 7));
  faces.emplace_back(Eigen::Vector3i(3, 7, 0));
  faces.emplace_back(Eigen::Vector3i(0, 7, 4));
  faces.emplace_back(Eigen::Vector3i(6, 5, 4));
  faces.emplace_back(Eigen::Vector3i(7, 6, 4));

  // Loading it.
  meshoui::Mesh mesh(vertices, faces);

  // Or reading from file
//  meshoui::Mesh mesh("../../data/Sphere.obj");

#ifdef MESHOUI_USE_VTK

  // VTKMesh.
  VTKMesh vtkmesh(mesh);

  // Building the different fields
  std::vector<double> faces_areas;
  faces_areas.reserve(mesh.n_faces());
  std::vector<double> face_centroid_z;
  face_centroid_z.reserve(mesh.n_faces());

  auto f_iter = mesh.faces_begin();
  for (; f_iter != mesh.faces_end(); ++f_iter) {
    faces_areas.push_back(mesh.data(*f_iter).GetSurfaceIntegral(meshoui::POLY_1));
    face_centroid_z.push_back(mesh.calc_face_centroid(*f_iter)[2]);
  }

  std::vector<double> vertices_z;
  vertices_z.reserve(mesh.n_vertices());
  auto v_iter = mesh.vertices_begin();
  for (; v_iter != mesh.vertices_end(); v_iter++) {
    vertices_z.push_back(mesh.point(*v_iter)[2]);
  }

  vtkmesh.AddField("faces_areas", faces_areas, VTKMesh::CELL);
  vtkmesh.AddField("face_centroid_z", face_centroid_z, VTKMesh::CELL);
  vtkmesh.AddField("vertices_z", vertices_z, VTKMesh::VERTEX);

  // Writing.
  vtkmesh.Write("mesh_field.vtp");

  // Visualization.
  vtkmesh.Visualize();

#endif

  return 0;
}

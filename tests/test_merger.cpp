// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "meshoui/meshoui.h"
#include "gtest/gtest.h"

using namespace meshoui;

// Test for checking the mesh merger.

TEST(meshoui_tests, merger) {

  #ifdef USE_VTK

  std::vector<Mesh*> meshes_two;
  meshes_two.reserve(2);
  std::vector<Mesh*> meshes_three;
  meshes_three.reserve(3);

  // First mesh.
  auto mesh_1 = new Mesh("../../data/Sphere.obj");
  meshes_two.push_back(mesh_1);
  meshes_three.push_back(mesh_1);
//  Show(mesh_1);

  // Second mesh.
//  Mesh mesh_2("../../data/Sphere.obj");
  auto mesh_2 = new Mesh("../../data/Sphere.obj");
  mesh_2->Translate(Vector3d(3., 0., 0.)); // Translation for seeing the two meshes.
  meshes_two.push_back(mesh_2);
  meshes_three.push_back(mesh_2);
//  Show(mesh_2);

  // Third mesh.
//  Mesh mesh_3("../../data/Sphere.obj");
  auto mesh_3 = new Mesh("../../data/Sphere.obj");
  mesh_3->Translate(Vector3d(0., 3., 0.)); // Translation for seeing the two meshes.
  meshes_three.push_back(mesh_3);
//  Show(mesh_3);

  // Merge - Two meshes.
  auto mesh_merged_two = Merge(meshes_two);
  assert(mesh_merged_two.n_faces() == (mesh_1->n_faces() + mesh_2->n_faces()));
  Show(mesh_merged_two);

  // Merge.
  auto mesh_merged_three = Merge(meshes_three);
  assert(mesh_merged_three.n_faces() == (mesh_1->n_faces() + mesh_2->n_faces() + mesh_3->n_faces()));
  Show(mesh_merged_three);

  auto merged = *mesh_1 + *mesh_2;
  assert(merged.n_faces() == (mesh_1->n_faces() + mesh_2->n_faces()));
  Show(merged);

  Write_VTK(mesh_merged_three, "Merged_mesh.vtp");

  #endif

}

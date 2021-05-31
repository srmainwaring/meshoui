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

  std::vector<Mesh*> meshes;
  meshes.reserve(3);

  // First mesh.
  auto mesh_1 = new Mesh("../../data/Sphere.obj");
  meshes.push_back(mesh_1);
//  Show(mesh_1);

  // Second mesh.
//  Mesh mesh_2("../../data/Sphere.obj");
  auto mesh_2 = new Mesh("../../data/Sphere.obj");
  mesh_2->Translate(Vector3d(3., 0., 0.)); // Translation for seeing the two meshes.
  meshes.push_back(mesh_2);
//  Show(mesh_2);

  // Third mesh.
//  Mesh mesh_3("../../data/Sphere.obj");
  auto mesh_3 = new Mesh("../../data/Sphere.obj");
  mesh_3->Translate(Vector3d(0., 3., 0.)); // Translation for seeing the two meshes.
  meshes.push_back(mesh_3);
//  Show(mesh_3);

  // Merge.
  auto mesh_merged = Merge(meshes);
  Show(mesh_merged);

  Write_VTK(mesh_merged, "Merged_mesh.vtp");

  #endif

}

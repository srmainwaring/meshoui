// ==========================================================================
// Helios
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for checking the symmetry of the mesh by a horizontal plane z = h.

int main() {

  // TODO: It is necessary to check this test. The feature of symmetrysing is not guaranteed.

  // Mesh.
  meshoui::Mesh mesh("../../data/Sphere.obj");

  // VTKMesh.
  VTKMesh vtkmesh(mesh);
  vtkmesh.AddFaceNormalField(&mesh);
  vtkmesh.AddVertexNormalField(&mesh);

  // Writing.
  vtkmesh.Write("Mesh.vtp");

  // Symmetry plane of equation z = 0;
  double height = 0;
  mesh.SymmetryHorizontalPlane(height);

  // Update properties.
  mesh.UpdateAllProperties();

  // Flipping both the face and vertex normals.
  mesh.FlipFaceNormals();
  mesh.FlipVertexNormals();

  // VTKMesh.
  VTKMesh vtkmesh_symmetrized(mesh);
  vtkmesh_symmetrized.AddFaceNormalField(&mesh);
  vtkmesh_symmetrized.AddVertexNormalField(&mesh);

  // Writing.
  vtkmesh_symmetrized.Write("Mesh_symmetrized.vtp");

  return 0;

}
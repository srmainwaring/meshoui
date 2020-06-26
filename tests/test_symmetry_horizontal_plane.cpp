
#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for checking the symmetry of the mesh by a horizontal plane z = h.

int main() {

  // Mesh.
  meshoui::Mesh mesh;
  mesh.Load("../../../Helios/docs/input_files/Sphere.obj");

#ifdef MESHOUI_USE_VTK
  // VTKMesh.
  VTKMesh vtkmesh(mesh);
  vtkmesh.AddFaceNormalField(&mesh);
  vtkmesh.AddVertexNormalField(&mesh);

  // Writing.
  vtkmesh.Write("Mesh.vtp");
#endif

  // Symmetry plane of equation z = 0;
  double height = 0;
  mesh.SymmetryHorizontalPlane(height);

  // Update properties.
  mesh.UpdateAllProperties();

  // Flipping both the face and vertex normals.
  mesh.FlipFaceNormals();
  mesh.FlipVertexNormals();

#ifdef MESHOUI_USE_VTK
  // VTKMesh.
  VTKMesh vtkmesh_symmetrized(mesh);
  vtkmesh_symmetrized.AddFaceNormalField(&mesh);
  vtkmesh_symmetrized.AddVertexNormalField(&mesh);

  // Writing.
  vtkmesh_symmetrized.Write("Mesh_symmetrized.vtp");
#endif

}

#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for checking the symmetry of the mesh by a horizontal plane z = h.

int main() {

  // Mesh.
  meshoui::Mesh mesh("../../../Helios/docs/input_files/Sphere.obj");

  // VTKMesh.
  VTKMesh vtkmesh(mesh);

  // Writing.
  vtkmesh.Write("Mesh.vtp");

  // Symmetry plane of equation z = 0;
  double height = 0;
  mesh.SymmetryHorizontalPlane(height);

  // VTKMesh.
  VTKMesh vtkmesh_symmetrized(mesh);

  // Writing.
  vtkmesh_symmetrized.Write("Mesh_symmetrized.vtp");

}
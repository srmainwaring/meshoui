
#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for loading a mesh input file.

int main() {

//    meshoui::Mesh mesh("../../docs/Input_files/Sphere_6200_faces.obj");
  meshoui::Mesh mesh("../../../Helios/docs/input_files/Sphere.obj");

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

    return 0;
}

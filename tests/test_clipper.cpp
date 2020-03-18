//
// Created by frongere on 18/03/2020.
//

#include "meshoui/meshoui.h"

using namespace meshoui;

int main () {

  Mesh mesh("data/Sphere.obj");

  Clipper<ClippingPlane> clipper;

  clipper.Apply(&mesh);




  VTKMesh vtk_mesh(mesh);

  vtk_mesh.Write("essai.vtp");


}

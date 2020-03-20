//
// Created by frongere on 18/03/2020.
//

#include "meshoui/meshoui.h"

using namespace meshoui;

int main () {

  Mesh mesh("data/Aker.obj");
//  Show(mesh);
  Write_VTK(mesh, "initial.vtp");

  // Remesh
  Remesher remesher;

  remesher.SetAngleDetectionThreshold(3);
  remesher.SetHausdorffParam(0.05);
  remesher.SetConstantEdgeSize(1.);

  remesher.Apply(&mesh, 1);
//  Show(mesh);
  Write_VTK(mesh, "remeshed_full.vtp");

  // Clip
  Clipper<ClippingPlane> clipper;
  clipper.Apply(&mesh);
//  Show(mesh);
  Write_VTK(mesh, "clipped.vtp");

  // Remesh after clipping
  remesher.Apply(&mesh, 1);
  remesher.Apply(&mesh, 1);
  Show(mesh);
  Write_VTK(mesh, "clipped_remeshed.vtp");

  Write_OBJ(mesh, "final.obj");

}

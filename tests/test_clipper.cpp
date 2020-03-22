//
// Created by frongere on 18/03/2020.
//

#include "meshoui/meshoui.h"

using namespace meshoui;

int main () {

  Mesh mesh("data/Ship.obj");
//  Show(mesh);
  Write_VTK(mesh, "initial.vtp");

  // Remesh
  Remesher remesher;

  remesher.SetAngleDetectionThreshold(5);
  remesher.SetHausdorffParam(0.08);
  remesher.SetConstantEdgeSize(2.5);

  remesher.Apply(&mesh);
//  Show(mesh);
  Write_VTK(mesh, "remeshed_full.vtp");

  // Clip
  Clipper<ClippingPlane> clipper;
  mesh = clipper.Apply(mesh);
//  Show(mesh);
  Write_VTK(mesh, "clipped.vtp");

  // Remesh after clipping

  remesher.SetConstantEdgeSize(5);
  remesher.Apply(&mesh);
  Show(mesh);
  Write_VTK(mesh, "clipped_remeshed.vtp");

  Write_OBJ(mesh, "final.obj");


//  clipper.ExtractClippedPolygonSet();




}

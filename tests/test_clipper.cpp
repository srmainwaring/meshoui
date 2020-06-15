// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for checking the mesh clipping and remeshing.

int main () {

  Mesh mesh;
  mesh.Load("../../data/Ship.obj");

//  Show(mesh);
  Write_VTK(mesh, "initial.vtp");

  // Remesh
  Remesher remesher;

  remesher.SetAngleDetectionThreshold(5);
  remesher.SetHausdorffParam(0.5);
  remesher.SetConstantEdgeSize(1);

  remesher.RemeshIt(mesh);
//  Show(mesh);
  Write_VTK(mesh, "remeshed_full.vtp");

  // Clip
  auto plane = std::make_shared<Plane>(Vector3d{0., 0., 0.}, Vector3d {0., 0., 1.}); // FIXME: pourquoi ce constructeur ne fonctionen pas ???
  Clipper<ClippingPlane> clipper(std::make_shared<ClippingPlane>(plane));
  clipper.ClipIt(mesh);
//  Show(mesh);
  Write_VTK(mesh, "clipped.vtp");

  // Remesh after clipping

  remesher.SetHausdorffParam(0.05);
  remesher.SetConstantEdgeSize(4);
  remesher.RemeshIt(mesh);
  remesher.RemeshIt(mesh);

  Show(mesh);
  Write_VTK(mesh, "clipped_remeshed.vtp");
  Write_OBJ(mesh, "final.obj");

//  clipper.ExtractClippedPolygonSet();

}

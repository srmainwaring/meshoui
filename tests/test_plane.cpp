//
// Created by frongere on 20/03/2020.
//


#include "meshoui/meshoui.h"

using namespace meshoui;

int main() {

  Vector3d origin(2, 0, 0);
//  Vector3d origin(0, 0, -2);


  Vector3d normal(1, 0, 1);
  normal.normalize();


  Plane plane(origin, normal);

  Vector3d proj = plane.GetPlaneFrameOrigin(); // {1, 0, 1}
  double alpha = GetAngleBetweenVectors({0., 0., 1.}, normal) * MU_180_PI; // 45


  // Transform validation
  Vector3d point(0, 0, 2);

  Vector3d p1 = plane.GetPointPositionInPlaneLocalFrame(point); // {-2*sqrt(2), 0, 0}
  Vector3d p2 = plane.GetPointPositionInPlaneLocalFrame(plane.GetPlaneFrameOrigin()); // {-sqrt(2), 0, 0}
  Vector3d p3 = plane.GetPointPositionInPlaneLocalFrame(plane.origin()); // {0, 0, 0}
  Vector3d p4 = plane.GetPointPositionInPlaneLocalFrame({2, 0, 2}); // {-sqrt(2), 0, sqrt(2)}

  // Features
  double offset = plane.GetPlaneOffset(); // sqrt(2)

  double d1 = plane.GetSignedDistanceToPoint({0., 0., 0.}); // -sqrt(2)
  double d2 = plane.GetSignedDistanceToPoint({2, 0, 0}); // 0
  double d3 = plane.GetSignedDistanceToPoint({0, 0, 2}); // 0
  double d4 = plane.GetSignedDistanceToPoint(plane.GetPlaneFrameOrigin()); // 0


  assert(plane.IsPointOnPlane({2, 0, 0}));
  assert(plane.IsPointOnPlane(origin));




  return 0;

}

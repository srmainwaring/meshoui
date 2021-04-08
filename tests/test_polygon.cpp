// ==========================================================================
// Helios
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "gtest/gtest.h"
#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for checking the use of polygons.

TEST(meshoui_tests, polygon) {

  Planar3DPolygon polygon(std::make_shared<Plane>());

  polygon.AddPoint({0, 0, 0});
  polygon.AddPoint({1, 0, 0});
  polygon.AddPoint({1, 1, 0});
  polygon.AddPoint({0, 1, 0});

  polygon.Update();

//  double area = polygon.GetArea();

}

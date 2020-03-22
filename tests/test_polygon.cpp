//
// Created by frongere on 21/03/2020.
//

#include "meshoui/meshoui.h"

using namespace meshoui;

int main() {

  Planar3DPolygon polygon(std::make_shared<Plane>());

  polygon.AddPoint({0, 0, 0});
  polygon.AddPoint({1, 0, 0});
  polygon.AddPoint({1, 1, 0});
  polygon.AddPoint({0, 1, 0});

  polygon.Update();

//  double area = polygon.GetArea();


  return 0;
}

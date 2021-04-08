// ==========================================================================
// Helios
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================
#include "gtest/gtest.h"
#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for checking the numerical surface integration.

// Definition of a function to integration.
class IntegrandTest : public IntegrandOnFace<double> {
 public:
  double Evaluate(const Vector3d &x) const override {
    return (2 - x(0) - 2 * x(1)); // f(x,y,z) = 2 - x - 2y.
  }
};

TEST(meshoui_tests, NumericalIntegration) {

  // Mesh.
  meshoui::Mesh mesh("../../data/Face_test_integration.obj");

  // Order of the integration.
  int order = 2;

  // Integrator.
  IntegrandTest myFunction;
  auto myIntegrator = integration<double>(&myFunction, order, &mesh);

  // Analytical integration.
  double analytical_result = 1. / 3.;

  // Numerical integration.
  auto f_iter = mesh.faces_begin();
  for (; f_iter != mesh.faces_end(); ++f_iter) {

    double numerical_result = myIntegrator.Compute(*f_iter);

    std::cout << "" << std::endl;
    std::cout << "Analytical result: " << analytical_result << std::endl;
    std::cout << "Numerical result: " << numerical_result << std::endl;
    std::cout << "Relative error: " << (100 * (numerical_result - analytical_result) / analytical_result) << std::endl;

    ASSERT_NEAR(numerical_result, analytical_result, 1e-12);

  }

}

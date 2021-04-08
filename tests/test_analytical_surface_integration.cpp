// ==========================================================================
// Helios
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "meshoui/meshoui.h"
#include "gtest/gtest.h"


using namespace meshoui;

// Test for checking the analytical surface integrations.

// All integral computations in a single integrand.
class Integrand : public IntegrandOnFace<VectorN> {
 public:
  Integrand() {}
  VectorN Evaluate(const Vector3d &x) const override {
    VectorN integrand = VectorN(16);
    integrand[0] = 1.;
    integrand[1] = x[0];
    integrand[2] = x[1];
    integrand[3] = x[2];
    integrand[4] = x[0] * x[0];
    integrand[5] = x[1] * x[1];
    integrand[6] = x[2] * x[2];
    integrand[7] = x[0] * x[1];
    integrand[8] = x[0] * x[2];
    integrand[9] = x[1] * x[2];
    integrand[10] = x[0] * x[0] * x[0];
    integrand[11] = x[1] * x[1] * x[1];
    integrand[12] = x[2] * x[2] * x[2];
    integrand[13] = x[0] * x[0] * x[1];
    integrand[14] = x[1] * x[1] * x[2];
    integrand[15] = x[2] * x[2] * x[0];
    return integrand;
  }
};


TEST(meshoui_tests, AnalitycalIntegration) {


  // Mesh.
  meshoui::Mesh mesh("../../data/Sphere.obj");

  // Number of integrations.
  int num_int = 16;

  // Order of the integration.
  int order = 5;

  // Error initialization.
  std::vector<std::vector<double>> error;
  error.reserve(num_int);
  std::vector<double> error_tmp;
  error_tmp.reserve(mesh.n_faces());
  for (int i = 0; i < num_int; ++i){
    error.push_back(error_tmp);
  }

  auto f_iter = mesh.faces_begin();
  for(; f_iter != mesh.faces_end(); ++f_iter) {

    // Numerical integrations.
    Integrand my_Integrand;
    auto my_Integrator = integration<VectorN>(&my_Integrand, order, &mesh);
    VectorN int_num = my_Integrator.Compute(*f_iter);

    // Analytical computations.
    VectorN int_ana = VectorN::Zero(num_int);
    int_ana[0] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_1);
    int_ana[1] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_X);
    int_ana[2] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Y);
    int_ana[3] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Z);
    int_ana[4] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_X2);
    int_ana[5] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Y2);
    int_ana[6] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Z2);
    int_ana[7] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_XY);
    int_ana[8] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_XZ);
    int_ana[9] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_YZ);
    int_ana[10] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_X3);
    int_ana[11] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Y3);
    int_ana[12] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Z3);
    int_ana[13] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_X2Y);
    int_ana[14] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Y2Z);
    int_ana[15] = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Z2X);

    // Absolute errors.
    for (int i = 0; i < num_int; ++i){
      error[i].push_back(abs(int_num[i] - int_ana[i]) / int_ana[i]);
      EXPECT_NEAR(int_num[i], int_ana[i], 1e-15);
    }
  }

  #ifdef USE_VTK
  // VTKMesh.
  VTKMesh vtkmesh(mesh);
  vtkmesh.AddField("01_POLY_1_error", error[0], VTKMesh::CELL);
  vtkmesh.AddField("02_POLY_X_error", error[1], VTKMesh::CELL);
  vtkmesh.AddField("03_POLY_Y_error", error[2], VTKMesh::CELL);
  vtkmesh.AddField("04_POLY_Z_error", error[3], VTKMesh::CELL);
  vtkmesh.AddField("05_POLY_X2_error", error[4], VTKMesh::CELL);
  vtkmesh.AddField("06_POLY_Y2_error", error[5], VTKMesh::CELL);
  vtkmesh.AddField("07_POLY_Z2_error", error[6], VTKMesh::CELL);
  vtkmesh.AddField("08_POLY_XY_error", error[7], VTKMesh::CELL);
  vtkmesh.AddField("09_POLY_XZ_error", error[8], VTKMesh::CELL);
  vtkmesh.AddField("10_POLY_YZ_error", error[9], VTKMesh::CELL);
  vtkmesh.AddField("11_POLY_X3_error", error[10], VTKMesh::CELL);
  vtkmesh.AddField("12_POLY_Y3_error", error[11], VTKMesh::CELL);
  vtkmesh.AddField("13_POLY_Z3_error", error[12], VTKMesh::CELL);
  vtkmesh.AddField("14_POLY_X2Y_error", error[13], VTKMesh::CELL);
  vtkmesh.AddField("15_POLY_Y2Z_error", error[14], VTKMesh::CELL);
  vtkmesh.AddField("16_POLY_Z2X_error", error[15], VTKMesh::CELL);

  // Writing.
  vtkmesh.Write("Comparison_Surface_integration_triangle.vtp");
  #endif

}

// ==========================================================================
// Helios
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "meshoui/meshoui.h"

using namespace meshoui;

// Test for checking the analytical surface integrations.

// POLY_1.
class Integrand_POLY_1 : public IntegrandOnFace<double> {
 public:
  Integrand_POLY_1() {}
  double Evaluate(const Vector3d &x) const override {
    return 1;
  }
};

// POLY_X, POLY_Y and POLY_Z.
class Integrand_POLY_LIN : public IntegrandOnFace<Vector3d> {
 public:
  Integrand_POLY_LIN() {}
  Vector3d Evaluate(const Vector3d &x) const override {
    return x;
  }
};

// POLY_XY, POLY_XZ, POLY_YZ, POLY_X2, POLY_Y2, POLY_Z2.
class Integrand_POLY_QUAD : public IntegrandOnFace<Matrix33> {
 public:
  Integrand_POLY_QUAD() {}
  Matrix33 Evaluate(const Vector3d &x) const override {
    return (x * x.transpose());
  }
};

int main() {

  //TODO: Verifier les integrales cubiques.

  // Mesh.
  meshoui::Mesh mesh("../../data/Sphere.obj");

  // Order of the integration.
  int order = 5;

  // Fields.
  std::vector<double> poly_1_err;
  poly_1_err.reserve(mesh.n_faces());
  std::vector<double> poly_X_err;
  poly_X_err.reserve(mesh.n_faces());
  std::vector<double> poly_Y_err;
  poly_Y_err.reserve(mesh.n_faces());
  std::vector<double> poly_Z_err;
  poly_Z_err.reserve(mesh.n_faces());
  std::vector<double> poly_XY_err;
  poly_XY_err.reserve(mesh.n_faces());
  std::vector<double> poly_XZ_err;
  poly_XZ_err.reserve(mesh.n_faces());
  std::vector<double> poly_YZ_err;
  poly_YZ_err.reserve(mesh.n_faces());
  std::vector<double> poly_X2_err;
  poly_X2_err.reserve(mesh.n_faces());
  std::vector<double> poly_Y2_err;
  poly_Y2_err.reserve(mesh.n_faces());
  std::vector<double> poly_Z2_err;
  poly_Z2_err.reserve(mesh.n_faces());

  // Numerical integration.
  auto f_iter = mesh.faces_begin();
  for(; f_iter != mesh.faces_end(); ++f_iter) {

    // POLY_1.
    Integrand_POLY_1 my_Integrand_POLY_1;
    auto my_Integrator_POLY_1 = integration<double>(&my_Integrand_POLY_1, order, &mesh);
    double poly_1_num = my_Integrator_POLY_1.Compute(*f_iter);
    double poly_1_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_1);
    poly_1_err.push_back(abs(poly_1_num - poly_1_ana));

    // POLY_X, POLY_Y and POLY_Z.
    Integrand_POLY_LIN my_Integrand_POLY_LIN;
    auto my_Integrator_POLY_LIN = integration<Vector3d>(&my_Integrand_POLY_LIN, order, &mesh);
    Vector3d poly_LIN_num = my_Integrator_POLY_LIN.Compute(*f_iter);
    double poly_X_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_X);
    double poly_Y_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Y);
    double poly_Z_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Z);
    poly_X_err.push_back(abs(poly_LIN_num[0] - poly_X_ana));
    poly_Y_err.push_back(abs(poly_LIN_num[1] - poly_Y_ana));
    poly_Z_err.push_back(abs(poly_LIN_num[2] - poly_Z_ana));

    // POLY_XY, POLY_XZ, POLY_YZ, POLY_X2, POLY_Y2, POLY_Z2.
    Integrand_POLY_QUAD my_Integrand_POLY_QUAD;
    auto my_Integrator_POLY_QUAD = integration<Matrix33>(&my_Integrand_POLY_QUAD, order, &mesh);
    Matrix33 poly_QUAD_num = my_Integrator_POLY_QUAD.Compute(*f_iter);
    double poly_XY_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_XY);
    double poly_XZ_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_XZ);
    double poly_YZ_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_YZ);
    double poly_X2_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_X2);
    double poly_Y2_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Y2);
    double poly_Z2_ana = mesh.data(*f_iter).GetSurfaceIntegral(POLY_Z2);
    poly_XY_err.push_back(abs(poly_QUAD_num(0, 1) - poly_XY_ana));
    poly_XZ_err.push_back(abs(poly_QUAD_num(0, 2) - poly_XZ_ana));
    poly_YZ_err.push_back(abs(poly_QUAD_num(1, 2) - poly_YZ_ana));
    poly_X2_err.push_back(abs(poly_QUAD_num(0, 0) - poly_X2_ana));
    poly_Y2_err.push_back(abs(poly_QUAD_num(1, 1) - poly_Y2_ana));
    poly_Z2_err.push_back(abs(poly_QUAD_num(2, 2) - poly_Z2_ana));
  }

  // VTKMesh.
  VTKMesh vtkmesh(mesh);
  vtkmesh.AddField("POLY_1_error", poly_1_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_X_error", poly_X_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_Y_error", poly_Y_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_Z_error", poly_Z_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_XY_error", poly_XY_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_XZ_error", poly_XZ_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_YZ_error", poly_YZ_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_X2_error", poly_X2_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_Y2_error", poly_Y2_err, VTKMesh::CELL);
  vtkmesh.AddField("POLY_Z2_error", poly_Z2_err, VTKMesh::CELL);

  // Writing.
  vtkmesh.Write("Comparison_Surface_integration_triangle.vtp");

  return 0;
}
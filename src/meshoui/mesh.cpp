// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include <iostream>
#include <unsupported/Eigen/EulerAngles>

#include "mesh.h"

namespace meshoui {

  Mesh::Mesh(const std::string &meshfile) {

    // Constructor of the class.

    Load(meshfile);
  }

  void Mesh::Load(const std::string &meshfile) {

    // This function loads the mesh file.

    if (!OpenMesh::IO::read_mesh(*this, meshfile)) {
      std::cerr << "Meshfile " << meshfile << " could not be read\n";
      exit(EXIT_FAILURE);
    }

    UpdateAllProperties();
  }

  void Mesh::UpdateAllProperties() {
    // Computation of the normal vectors of faces, vertices and half edges.
    update_normals();

    for (FaceIter f_iter = faces_begin(); f_iter != faces_end(); ++f_iter) {

      // Face centroids.
      data(*f_iter).SetCenter(calc_face_centroid(*f_iter));

      // Face areas.
      CalcFacePolynomialIntegrals(*f_iter);
    }

  }

  void Mesh::CalcFacePolynomialIntegrals(const Mesh::FaceHandle &fh) {

    using Array = Eigen::Array<double, 3, 1>;

    Array P0, P1, P2;
    Array t0, t1, t2;
    Array f1, f2, f3;
    Array g0, g1, g2;
    Point e1, e2, cp;

    double delta;

    // Getting one half-edge handle of the current face
    auto heh = halfedge_handle(fh);

    // Getting the origin vertex of heh
    P0 = point(from_vertex_handle(heh));

    heh = next_halfedge_handle(heh);
    P1 = point(from_vertex_handle(heh));

    heh = next_halfedge_handle(heh);
    P2 = point(from_vertex_handle(heh));

    e1 = P1 - P0;
    e2 = P2 - P0;
    cp = cross(e1, e2);
    delta = cp.norm();

    // Factorized terms.
    t0 = P0 + P1;
    t1 = P0 * P0;
    t2 = t1 + P1 * t0;

    f1 = t0 + P2;
    f2 = t2 + P2 * f1;
    f3 = P0 * t1 + P1 * t2 + P2 * f2;

    g0 = f2 + P0 * (f1 + P0);
    g1 = f2 + P1 * (f1 + P1);
    g2 = f2 + P2 * (f1 + P2);

    // Extended Eberly's Formulas.
    // Surface integrals are transformed into contour integrals.
    data(fh).SetSurfaceIntegral(POLY_1, delta / 2.);

    data(fh).SetSurfaceIntegral(POLY_X, delta * f1[0] / 6.);
    data(fh).SetSurfaceIntegral(POLY_Y, delta * f1[1] / 6.);
    data(fh).SetSurfaceIntegral(POLY_Z, delta * f1[2] / 6.);

    data(fh).SetSurfaceIntegral(POLY_XY, delta * (f1[0] * f1[1] + P0[0] * P0[1] + P1[0] * P1[1] + P2[0] * P2[1]) / 24.);
    data(fh).SetSurfaceIntegral(POLY_XZ, delta * (f1[0] * f1[2] + P0[0] * P0[2] + P1[0] * P1[2] + P2[0] * P2[2]) / 24.);
    data(fh).SetSurfaceIntegral(POLY_YZ, delta * (f1[1] * f1[2] + P0[1] * P0[2] + P1[1] * P1[2] + P2[1] * P2[2]) / 24.);

    data(fh).SetSurfaceIntegral(POLY_X2, delta * f2[0] / 12.);
    data(fh).SetSurfaceIntegral(POLY_Y2, delta * f2[1] / 12.);
    data(fh).SetSurfaceIntegral(POLY_Z2, delta * f2[2] / 12.);

    data(fh).SetSurfaceIntegral(POLY_X3, delta * f3[0] / 20.);
    data(fh).SetSurfaceIntegral(POLY_Y3, delta * f3[1] / 20.);
    data(fh).SetSurfaceIntegral(POLY_Z3, delta * f3[2] / 20.);

    data(fh).SetSurfaceIntegral(POLY_X2Y, delta * (P0[1] * g0[0] + P1[1] * g1[0] + P2[1] * g2[0]) / 60.);
    data(fh).SetSurfaceIntegral(POLY_Y2Z, delta * (P0[2] * g0[1] + P1[2] * g1[1] + P2[2] * g2[1]) / 60.);
    data(fh).SetSurfaceIntegral(POLY_Z2X, delta * (P0[0] * g0[2] + P1[0] * g1[2] + P2[0] * g2[2]) / 60.);

  }

  void Mesh::SymmetryHorizontalPlane(const double &height) {
    // Symmetry of the mesh.
    for (VertexIter v_iter = vertices_begin(); v_iter != vertices_end(); ++v_iter) {
      point(*v_iter)[2] = 2 * height - point(*v_iter)[2];
    }
  }

  void Mesh::Translate(const Vector3d &t) {
    for (VertexIter v_iter = vertices_begin(); v_iter != vertices_end(); ++v_iter) {
      point(*v_iter) += t;
    }
    UpdateAllProperties();
  }

  void Mesh::Rotate(double phi, double theta, double psi) {
    auto euler_angles = Eigen::EulerAngles<double, Eigen::EulerSystemZYX>(psi, theta,
                                                                          phi); // FIXME: verifier qu'on a le resultat escompte

    auto vh_iter = vertices_begin();
    for (; vh_iter != vertices_end(); ++vh_iter) {
      point(*vh_iter) = euler_angles * point(*vh_iter);
    }

    UpdateAllProperties();
  }

  void Mesh::FlipFaceNormals() {
    for (FaceIter f_iter = faces_begin(); f_iter != faces_end(); ++f_iter) {
      Vector3d face_normal = normal(*f_iter);
      set_normal(*f_iter, -face_normal);
    }
  }

  void Mesh::FlipVertexNormals() {
    // FIXME: ok, si on symmetrize les normales doivent etre flippees mais ne peut on pas plutot recalculer en appelant
    // un update ?
    for (VertexIter v_iter = vertices_begin(); v_iter != vertices_end(); ++v_iter) {
      Vector3d vertex_normal = normal(*v_iter);
      set_normal(*v_iter, -vertex_normal);
    }
  }

  void Mesh::FlipHalfedgeNormals() {
    for (HalfedgeIter he_iter = halfedges_begin(); he_iter != halfedges_end(); ++he_iter) {
      Vector3d halfedge_normal = normal(*he_iter);
      set_normal(*he_iter, -halfedge_normal);
    }
  }

  void Mesh::FlipAllNormals() {
    FlipFaceNormals();
    FlipVertexNormals();
    FlipHalfedgeNormals();
  }

} // end namespace meshoui

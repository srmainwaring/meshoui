// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include <iostream>

#include "mesh.h"
#include "math.h"

namespace meshoui {

//  Mesh::Mesh(const std::string &meshfile) {
//
//    // Constructor of the class.
//
//    Load(meshfile);
//  }

  void Mesh::Load(const std::vector<Vector3d> &vertices, const std::vector<Eigen::VectorXi> &faces) {
    //load vertices
    VertexHandle vertexHandle[vertices.size()];
    int iv = 0;
    for (const auto &vertex : vertices) {
      vertexHandle[iv] = add_vertex(vertex);
      iv++;
    }
    //load (triangular) faces
    for (const auto& face : faces) {
      std::vector<VertexHandle> face_vhandles;
      face_vhandles.clear();
      face_vhandles.push_back(vertexHandle[face[0]]);
      face_vhandles.push_back(vertexHandle[face[1]]);
      face_vhandles.push_back(vertexHandle[face[2]]);
      add_face(face_vhandles);
    }

      UpdateAllProperties();
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

    // This function updates some properties of faces and vertices (normals, centroids, face areas).

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

    // This function computes the polynomial surface integrals over the faces.

    // Getting one half-edge handle of the current face
    auto heh = halfedge_handle(fh);

    // Getting the origin vertex of heh
    Vector3d P0 = point(from_vertex_handle(heh));

    heh = next_halfedge_handle(heh);
    Vector3d P1 = point(from_vertex_handle(heh));

    heh = next_halfedge_handle(heh);
    Vector3d P2 = point(from_vertex_handle(heh));

    Vector3d e1 = P1 - P0;
    Vector3d e2 = P2 - P0;
    Vector3d cp = cross(e1, e2);
    double delta = cp.norm();

    // My Extended Eberly's Formulas.
    // Surface integrals are transformed into contour integrals.
    data(fh).SetSurfaceIntegral(AREA, delta / 2.);

  }

  void Mesh::SymmetryHorizontalPlane(const double &height) {

    // This function applies a symmetry by a plane of equation z = h.

    // Symmetry of the mesh.
    for (VertexIter v_iter = vertices_begin(); v_iter != vertices_end(); ++v_iter) {
      point(*v_iter)[2] = 2 * height - point(*v_iter)[2];
    }

  }

  void Mesh::FlipFaceNormals() {

    // This function flips the face normals.

    for (FaceIter f_iter = faces_begin(); f_iter != faces_end(); ++f_iter) {
      Vector3d face_normal = normal(*f_iter);
      set_normal(*f_iter, -face_normal);
    }

  }

  void Mesh::FlipVertexNormals() {

    // This function flips the vertex normals.

    for (VertexIter v_iter = vertices_begin(); v_iter != vertices_end(); ++v_iter) {
      Vector3d vertex_normal = normal(*v_iter);
      set_normal(*v_iter, -vertex_normal);
    }

  }

  void Mesh::FlipHalfedgeNormals() {

    // This function flips the halfegde normals.

    for (HalfedgeIter he_iter = halfedges_begin(); he_iter != halfedges_end(); ++he_iter) {
      Vector3d halfedge_normal = normal(*he_iter);
      set_normal(*he_iter, -halfedge_normal);
    }

  }

  void Mesh::FlipAllNormals() {

    // This function flips the normals (face, vertice, halfedge).

    FlipFaceNormals();
    FlipVertexNormals();
    FlipHalfedgeNormals();
  }

} // end namespace meshoui

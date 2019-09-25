// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include <iostream>

#include "mesh.h"

namespace meshoui {

    mesh::mesh(std::string meshfile) {

        // Constructor of the class.

        Load(std::move(meshfile));
    }

    void mesh::Load(std::string meshfile) {

        // This function loads the mesh file.

        if (!OpenMesh::IO::read_mesh(*this, meshfile)) {
            std::cerr << "Meshfile " << meshfile << " could not be read\n";
            exit(1);
        }

        UpdateAllProperties();
    }

    void mesh::UpdateAllProperties() {

        // This function updates all properties of faces and vertices (normals, centroids, surface integrals).

        // Computation of normal vectors and centroids.
        UpdateBaseProperties();

        // Computation of surface polynomial integrals.
        UpdateFacesPolynomialIntegrals();

    }

    void mesh::UpdateBaseProperties() {

        // This function computes the normal vectors everywhere and the centroid of faces.

        // Computation of the normal vectors of faces, vertices and half edges.
        update_normals(); // Update normals for both faces and vertices

        // Update face's center properties.
        Point center;
        for (FaceIter f_iter = faces_begin(); f_iter != faces_end(); ++f_iter) {
            data(*f_iter).SetCenter(calc_face_centroid(*f_iter));
        }

    }

    void mesh::UpdateFacesPolynomialIntegrals() {

        // This function updates the computations of the polynomial surface integrals.

        for (FaceIter f_iter = faces_begin(); f_iter != faces_end(); ++f_iter) {
            CalcFacePolynomialIntegrals(*f_iter);
        }
    }

    void mesh::CalcFacePolynomialIntegrals(const mesh::FaceHandle &fh) {

        // This function computes the polynomial surface integrals over the faces.

        typedef OpenMesh::Vec3d Point;

        Point P0, P1, P2;
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

        // My Extended Eberly's Formulas.
        // Surface integrals are transformed into contour integrals.
        data(fh).SetSurfaceIntegral(POLY_1, delta / 2.);

    }



} // end namespace meshoui

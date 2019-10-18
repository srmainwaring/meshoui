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

        Vector3d P0, P1, P2;
        Vector3d e1, e2, cp;

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
        data(fh).SetSurfaceIntegral(AREA, delta / 2.);

    }

} // end namespace meshoui

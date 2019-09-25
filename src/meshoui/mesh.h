// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_LOADER_H
#define MESHOUI_LOADER_H

#include "meshTraits.h"

namespace meshoui {

    class mesh : public OpenMesh::TriMesh_ArrayKernelT<meshTraits> { // mesh must be a triangular mesh.

    public:

        /// Constructor of the class.
        mesh() = default;

        /// Constructor of the class.
        explicit mesh(std::string meshfile);

        /// This function loads the mesh file.
        void Load(std::string meshfile);

        /// This function updates all properties of faces and vertices (normals, centroids, surface integrals).
        void UpdateAllProperties();

    private:

        /// This function computes the normal vectors everywhere and the centroid of faces.
        void UpdateBaseProperties();

        /// This function updates the computations of the polynomial surface integrals.
        void UpdateFacesPolynomialIntegrals();

        /// Computes triangular faces surface integration of some polynomial integrands using analytical formulas
        /// established by transforming surface integrals into contour integrals and deriving analytical expressions.
        /// Extended from Eberly... https://d-ice.gitlab.host/common/technical_reports/mesh-integrals
        void CalcFacePolynomialIntegrals(const mesh::FaceHandle &fh);

    };

} // end namespace meshoui

#endif // MESHOUI_LOADER_H

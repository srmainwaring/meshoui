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

    /**
    * Class for dealing with OpenMesh structures.
    */
    class mesh : public OpenMesh::TriMesh_ArrayKernelT<meshTraits> { // Mesh must be a triangular mesh.

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

        /// Convert an OpenMesh point into a vector.
        template <class Vector>
        inline Vector OpenMeshPointToVector3d(const mesh::Point &point) {
            return Vector(point[0], point[1], point[2]); // Always gives a FRyDoM vector expressed in NWU
        }

    };

} // end namespace meshoui

#endif // MESHOUI_LOADER_H
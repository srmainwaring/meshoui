// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_MESHTRAITS_H
#define MESHOUI_MESHTRAITS_H

#include <iostream>

#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>

namespace meshoui {

    enum IntegrandType {
        // TODO: voir http://www.drdobbs.com/when-enum-just-isnt-enough-enumeration-c/184403955 pour une meilleure
        // gestion des enums
        UNDEFINED_INTEGRAND,
        AREA
    };

struct meshTraits : public OpenMesh::DefaultTraits {
        typedef OpenMesh::Vec3d Point;

        VertexAttributes(OpenMesh::Attributes::Normal |
                                 OpenMesh::Attributes::Status
        );

        FaceAttributes(OpenMesh::Attributes::Normal |
                               OpenMesh::Attributes::TAGGED |
                               OpenMesh::Attributes::FEATURE |
                               OpenMesh::Attributes::Status
        );

        HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge |
                                   OpenMesh::Attributes::Status |
                                   OpenMesh::Attributes::TAGGED
        );

        EdgeAttributes(OpenMesh::Attributes::Status);

        EdgeTraits
        {
            private:
            double m_length = 0.0;

            public:
            EdgeT() {}

            const double GetLength() const {
                return m_length;
            }

            void SetLength(double length) {
                m_length = length;
            }

        };

        FaceTraits
        {
            private:
            Point m_center = {0.0, 0.0, 0.0};

            struct SurfaceIntegrals {
                double m_int_1 = 0.; // TODO: calculer int_1, pas fait encore
            };

            SurfaceIntegrals m_integrals;

            public:
            FaceT() {}

            const Point &Center() const { return m_center; }

            /// This function sets the position of the face centroid.
            void SetCenter(const Point &center) { m_center = center; }

            /// This function gives the surface integral of a face.
            const double GetSurfaceIntegral(IntegrandType type) const { // TODO: abandonner les enums pour les integrandes et preferer les accessors voir mieux, des fonctors...
                switch (type) {
                    case AREA:
                        return m_integrals.m_int_1; // This is the surface area...
                    case UNDEFINED_INTEGRAND:
                        std::cerr << "Cannot return value of an UNDEFINED_INTEGRAND" << std::endl;
                        break;
                    default:
                        throw std::runtime_error("Cannot return value of an UNDEFINED_INTEGRAND");
                }
            }

            void SetSurfaceIntegral(IntegrandType type, const double &val) {
                switch (type) {
                    case AREA:
                        m_integrals.m_int_1 = val;
                        break;
                    case UNDEFINED_INTEGRAND:
                        std::cerr << "Cannot return value of an UNDEFINED_INTEGRAND" << std::endl;
                        break;
                    default:
                        throw std::runtime_error("Cannot return value of an UNDEFINED_INTEGRAND");
                }
            }

        };

    };

} // end namespace meshoui

#endif //MESHOUI_MESHTRAITS_H

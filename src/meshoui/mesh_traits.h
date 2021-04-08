// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_MESH_TRAITS_H
#define MESHOUI_MESH_TRAITS_H

#include <iostream>

// FIXME: Ces include devraient plutot se trouver dans mesh.h !!
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>

namespace meshoui {

  enum IntegrandType {
    // TODO: voir http://www.drdobbs.com/when-enum-just-isnt-enough-enumeration-c/184403955 pour une meilleure
    // gestion des enums
    UNDEFINED_INTEGRAND,
    POLY_1,
    POLY_X,
    POLY_Y,
    POLY_Z,
    POLY_YZ,
    POLY_XZ,
    POLY_XY,
    POLY_X2,
    POLY_Y2,
    POLY_Z2,
    POLY_X3,
    POLY_Y3,
    POLY_Z3,
    POLY_X2Y,
    POLY_Y2Z,
    POLY_Z2X,
    POLY_Y2X
  };

  struct meshTraits : public OpenMesh::DefaultTraits {

    // Use of Eigen with OpenMesh.
    using Point = Eigen::Vector3d;
    using Normal = Eigen::Vector3d;
    using TexCoord2D = Eigen::Vector2d;

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
        double m_int_1 = 0.;

        double m_int_x = 0.;
        double m_int_y = 0.;
        double m_int_z = 0.;

        double m_int_yz = 0.;
        double m_int_xz = 0.;
        double m_int_xy = 0.;

        double m_int_x2 = 0.;
        double m_int_y2 = 0.;
        double m_int_z2 = 0.;

        double m_int_x3 = 0.;
        double m_int_y3 = 0.;
        double m_int_z3 = 0.;

        double m_int_x2y = 0.;
        double m_int_y2z = 0.;
        double m_int_z2x = 0.;
      };

      SurfaceIntegrals m_integrals;

     public:
      FaceT() {}

      const Point &Center() const { return m_center; }

      /// This function sets the position of the face centroid.
      void SetCenter(const Point &center) { m_center = center; }

      /// This function gives the surface integral of a face.
      const double GetSurfaceIntegral(IntegrandType type) const {
        // TODO: abandonner les enums pour les integrandes et preferer les accessors voir mieux, des fonctors...
        double si;
        switch (type) {
          case POLY_1:
            si = m_integrals.m_int_1; // This is the surface area...
            break;
          case POLY_X:
            si = m_integrals.m_int_x;
            break;
          case POLY_Y:
            si = m_integrals.m_int_y;
            break;
          case POLY_Z:
            si = m_integrals.m_int_z;
            break;
          case POLY_YZ:
            si = m_integrals.m_int_yz;
            break;
          case POLY_XZ:
            si = m_integrals.m_int_xz;
            break;
          case POLY_XY:
            si = m_integrals.m_int_xy;
            break;
          case POLY_X2:
            si = m_integrals.m_int_x2;
            break;
          case POLY_Y2:
            si = m_integrals.m_int_y2;
            break;
          case POLY_Z2:
            si = m_integrals.m_int_z2;
            break;
          case POLY_X3:
            si = m_integrals.m_int_x3;
            break;
          case POLY_Y3:
            si = m_integrals.m_int_y3;
            break;
          case POLY_Z3:
            si = m_integrals.m_int_z3;
            break;
          case POLY_X2Y:
            si = m_integrals.m_int_x2y;
            break;
          case POLY_Y2Z:
            si = m_integrals.m_int_y2z;
            break;
          case POLY_Z2X:
            si = m_integrals.m_int_z2x;
            break;
          case UNDEFINED_INTEGRAND:
            std::cerr << "Cannot return value of an UNDEFINED_INTEGRAND" << std::endl;
            break;
          default:
            throw std::runtime_error("Cannot return value of an UNDEFINED_INTEGRAND");
        }
        return si;
      }

      void SetSurfaceIntegral(IntegrandType type, const double &val) {
        switch (type) {
          case POLY_1:
            m_integrals.m_int_1 = val;
            break;
          case POLY_X:
            m_integrals.m_int_x = val;
            break;
          case POLY_Y:
            m_integrals.m_int_y = val;
            break;
          case POLY_Z:
            m_integrals.m_int_z = val;
            break;
          case POLY_YZ:
            m_integrals.m_int_yz = val;
            break;
          case POLY_XZ:
            m_integrals.m_int_xz = val;
            break;
          case POLY_XY:
            m_integrals.m_int_xy = val;
            break;
          case POLY_X2:
            m_integrals.m_int_x2 = val;
            break;
          case POLY_Y2:
            m_integrals.m_int_y2 = val;
            break;
          case POLY_Z2:
            m_integrals.m_int_z2 = val;
            break;
          case POLY_X3:
            m_integrals.m_int_x3 = val;
            break;
          case POLY_Y3:
            m_integrals.m_int_y3 = val;
            break;
          case POLY_Z3:
            m_integrals.m_int_z3 = val;
            break;
          case POLY_X2Y:
            m_integrals.m_int_x2y = val;
            break;
          case POLY_Y2Z:
            m_integrals.m_int_y2z = val;
            break;
          case POLY_Z2X:
            m_integrals.m_int_z2x = val;
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

#endif //MESHOUI_MESH_TRAITS_H

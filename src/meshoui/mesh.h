// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_LOADER_H
#define MESHOUI_LOADER_H

#include "meshTraits.h"
#include <OpenMesh/Core/Utils/PropertyManager.hh>

namespace meshoui {

    // OpenMesh handles.
    using FaceHandle = OpenMesh::FaceHandle;
    using VertexHandle = OpenMesh::VertexHandle;
    using HalfedgeHandle = OpenMesh::HalfedgeHandle;
    using EdgeHandle = OpenMesh::EdgeHandle;

    /**
    * Class for dealing with OpenMesh structures.
    */
    class Mesh : public OpenMesh::TriMesh_ArrayKernelT<meshTraits> { // Mesh must be a triangular mesh.

    public:

        /// Constructor of the class.
        Mesh() = default;

        /// Constructor of the class.
        explicit Mesh(std::string meshfile);

        /// This function loads the mesh file.
        void Load(std::string meshfile);

        /// This function updates all properties of faces and vertices (normals, centroids, surface integrals).
        void UpdateAllProperties();

        /// This function returns a property about the faces of the mesh.
        template<typename T>
        OpenMesh::PropertyManager<OpenMesh::FPropHandleT<T>, meshoui::Mesh> GetFaceProperty(const char* name){
          return OpenMesh::getProperty<meshoui::FaceHandle, T>(*this, name);
        }

        /// This function returns a property about the vertices of the mesh.
        template<typename T>
        OpenMesh::PropertyManager<OpenMesh::VPropHandleT<T>, meshoui::Mesh> GetVertexProperty(const char* name){
          return OpenMesh::getProperty<meshoui::VertexHandle, T>(*this, name);
        }

        /// This function returns a property about the edges of the mesh.
        template<typename T>
        OpenMesh::PropertyManager<OpenMesh::EPropHandleT<T>, meshoui::Mesh> GetEdgeProperty(const char* name){
          return OpenMesh::getProperty<meshoui::EdgeHandle, T>(*this, name);
        }

      /// This function creates a property about the faces of the mesh.
      template<typename T>
      void CreateFaceProperty(const char* name){
        OpenMesh::getOrMakeProperty<meshoui::FaceHandle, T>(*this, name);
      }

    private:

        /// This function computes the normal vectors everywhere and the centroid of faces.
        void UpdateBaseProperties();

        /// This function updates the computations of the polynomial surface integrals.
        void UpdateFacesPolynomialIntegrals();

        /// Computes triangular faces surface integration of some polynomial integrands using analytical formulas
        /// established by transforming surface integrals into contour integrals and deriving analytical expressions.
        /// Extended from Eberly... https://d-ice.gitlab.host/common/technical_reports/mesh-integrals
        void CalcFacePolynomialIntegrals(const Mesh::FaceHandle &fh);

    };

  /// Convert an OpenMesh point into a vector.
  template <class Vector>
  inline Vector OpenMeshPointToVector3d(const Mesh::Point &point) {
    return Vector(point[0], point[1], point[2]); // Always gives a FRyDoM vector expressed in NWU
  }

} // end namespace meshoui

#endif // MESHOUI_LOADER_H

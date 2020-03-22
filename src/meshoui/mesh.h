// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_LOADER_H
#define MESHOUI_LOADER_H

#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <unordered_map>

#include "maths.h"
#include "meshTraits.h"

#include "polygon.h"


namespace meshoui {

  // OpenMesh handles.
  using FaceHandle = OpenMesh::FaceHandle;
  using VertexHandle = OpenMesh::VertexHandle;
  using HalfedgeHandle = OpenMesh::HalfedgeHandle;
  using EdgeHandle = OpenMesh::EdgeHandle;
  using FPropHandleTDouble = OpenMesh::FPropHandleT<double>;
  using VPropHandleTDouble = OpenMesh::VPropHandleT<double>;
  using EPropHandleTDouble = OpenMesh::EPropHandleT<double>;
  using VertexIter = OpenMesh::PolyConnectivity::VertexIter;
  using FaceIter = OpenMesh::PolyConnectivity::FaceIter;

  template<class T>
  using FPropHandleTMap = OpenMesh::FPropHandleT<std::unordered_map<VertexHandle, std::array<T, 2>>>;

  /**
  * Class for dealing with OpenMesh structures.
  */
  class Mesh : public OpenMesh::TriMesh_ArrayKernelT<meshTraits> { // Mesh must be a triangular mesh.

   public:

//    using Transform = Eigen::Transform<double, 3, Eig>; // TODO: faire quelque chose de dedie dans mathutils

    // OpenMesh property managers.
    template<typename T>
    using FaceProperty = OpenMesh::PropertyManager<OpenMesh::FPropHandleT<T>, meshoui::Mesh>;
    template<typename T>
    using VertexProperty = OpenMesh::PropertyManager<OpenMesh::VPropHandleT<T>, meshoui::Mesh>;
    template<typename T>
    using EdgeProperty = OpenMesh::PropertyManager<OpenMesh::EPropHandleT<T>, meshoui::Mesh>;

    /// Default constructor.
    Mesh() = default;

    /// Constructor of the class.
    explicit Mesh(const std::string &meshfile);

    /// This function loads the mesh file.
    void Load(const std::string &meshfile);

    /// This function updates some properties of faces and vertices (normals, centroids, face areas).
    void UpdateAllProperties();

    /// This function returns a property about the faces of the mesh.
    template<typename T>
    FaceProperty<T> GetFaceProperty(const char *name) {
      return OpenMesh::getProperty<meshoui::FaceHandle, T>(*this, name);
    }

    /// This function returns a property about the vertices of the mesh.
    template<typename T>
    VertexProperty<T> GetVertexProperty(const char *name) {
      return OpenMesh::getProperty<meshoui::VertexHandle, T>(*this, name);
    }

    /// This function returns a property about the edges of the mesh.
    template<typename T>
    EdgeProperty<T> GetEdgeProperty(const char *name) {
      return OpenMesh::getProperty<meshoui::EdgeHandle, T>(*this, name);
    }

    /// This function creates a property about the faces of the mesh.
    template<typename T>
    FaceProperty<T> CreateFaceProperty(const char *name) {
      return OpenMesh::getOrMakeProperty<meshoui::FaceHandle, T>(*this, name);
    }

    /// This function creates a property about the vertices of the mesh.
    template<typename T>
    VertexProperty<T> CreateVertexProperty(const char *name) {
      return OpenMesh::getOrMakeProperty<meshoui::VertexHandle, T>(*this, name);
    }

    /// This function creates a property about the edges of the mesh.
    template<typename T>
    EdgeProperty<T> CreateEdgeProperty(const char *name) {
      return OpenMesh::getOrMakeProperty<meshoui::EdgeHandle, T>(*this, name);
    }

    /// This function creates a temporary property about the faces of the mesh.
    template<typename T>
    FaceProperty<T> CreateTemporaryFaceProperty() {
      return OpenMesh::makeTemporaryProperty<meshoui::FaceHandle, T>(*this);
    }

    /// This function creates a temporary property about the vertices of the mesh.
    template<typename T>
    VertexProperty<T> CreateTemporaryVertexProperty() {
      return OpenMesh::makeTemporaryProperty<meshoui::VertexHandle, T>(*this);
    }

    /// This function creates a temporary property about the edges of the mesh.
    template<typename T>
    EdgeProperty<T> CreateTemporaryEdgeProperty() {
      return OpenMesh::makeTemporaryProperty<meshoui::EdgeHandle, T>(*this);
    }

    /// This function applies a symmetry by a plane of equation z = h.
    void SymmetryHorizontalPlane(const double &height);

    /// This function translates the mesh.
    void Translate(const Vector3d& t);

    /// This function rotates the mesh, based on Cardan angles
    void Rotate(double phi, double theta, double psi);

    //TODO:
//    void AffineTransform(const Eigen::Transform<double, 3, Eigen::Affine>& transform);

    bool IsWatertight() const;

//    // FIXME: Pas forcement ici...
//    Planar3DPolygonSet ExtractIntersectionPolygons(std::shared_ptr<Plane> plane);



    /// This function flips the face normals.
    void FlipFaceNormals();

    /// This function flips the vertex normals.
    void FlipVertexNormals();

    /// This function flips the halfegde normals.
    void FlipHalfedgeNormals();

    /// This function flips the normals (face, vertice, halfedge).
    void FlipAllNormals();

   private:

    /// Computes triangular faces surface integration of some polynomial integrands using analytical formulas
    /// established by transforming surface integrals into contour integrals and deriving analytical expressions.
    /// Extended from Eberly... https://d-ice.gitlab.host/common/technical_reports/mesh-integrals
    void CalcFacePolynomialIntegrals(const Mesh::FaceHandle &fh);

  };

} // end namespace meshoui

#endif // MESHOUI_LOADER_H

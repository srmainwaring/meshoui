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
#include "mesh_traits.h"

#include "polygon.h"

namespace meshoui {

  /// OpenMesh handles.
  using FaceHandle = OpenMesh::FaceHandle;
  using VertexHandle = OpenMesh::VertexHandle;
  using HalfedgeHandle = OpenMesh::HalfedgeHandle;
  using EdgeHandle = OpenMesh::EdgeHandle;
  using FPropHandleTDouble = OpenMesh::FPropHandleT<double>;
  using VPropHandleTDouble = OpenMesh::VPropHandleT<double>;
  using EPropHandleTDouble = OpenMesh::EPropHandleT<double>;
  using VertexIter = OpenMesh::PolyConnectivity::VertexIter;
  using FaceIter = OpenMesh::PolyConnectivity::FaceIter;

  template<class PropertyType>
  using FPropHandleTMap = OpenMesh::FPropHandleT<std::unordered_map<VertexHandle, std::array<PropertyType, 2>>>;

  /**
  * Class for dealing with OpenMesh structures.
  */
  class Mesh : public OpenMesh::TriMesh_ArrayKernelT<meshTraits> { // Mesh must be a triangular mesh.

   public:

    /// OpenMesh property managers.
    template<typename PropertyType>
    using FaceProperty = OpenMesh::PropertyManager<OpenMesh::FPropHandleT<PropertyType>, Mesh>;
    template<typename PropertyType>
    using VertexProperty = OpenMesh::PropertyManager<OpenMesh::VPropHandleT<PropertyType>, Mesh>;
    template<typename PropertyType>
    using EdgeProperty = OpenMesh::PropertyManager<OpenMesh::EPropHandleT<PropertyType>, Mesh>;

    Mesh() = default;

    /// Constructor of the class.
    explicit Mesh(const std::string &meshfile);

    /// Constructor of the class.
    Mesh(const Mesh& mesh);

    /// Constructor of the class.
    Mesh(const std::vector<Vector3d> &vertices, const std::vector<Eigen::VectorXi> & faces);

    /// Constructor of the class.
    void Load(const std::string &meshfile);

    /// Constructor of the class.
    void Load(const std::vector<Vector3d> &vertices, const std::vector<Eigen::VectorXi> & faces);

    /// This method updates some properties of faces and vertices (normals, centroids, face areas).
    void UpdateAllProperties();

    /// This method returns a property about the faces of the mesh.
    template<typename PropertyType>
    FaceProperty<PropertyType> GetFaceProperty(const char* name){
      return OpenMesh::getProperty<meshoui::FaceHandle, PropertyType>(*this, name);
    }

    /// This method returns a property about the vertices of the mesh.
    template<typename PropertyType>
    VertexProperty<PropertyType> GetVertexProperty(const char *name) {
      return OpenMesh::getProperty<VertexHandle, PropertyType>(*this, name);
    }

    /// This method returns a property about the edges of the mesh.
    template<typename PropertyType>
    EdgeProperty<PropertyType> GetEdgeProperty(const char *name) {
      return OpenMesh::getProperty<EdgeHandle, PropertyType>(*this, name);
    }

    /// This method creates a property about the faces of the mesh.
    template<typename PropertyType>
    FaceProperty<PropertyType> CreateFaceProperty(const char *name) {
      return OpenMesh::getOrMakeProperty<FaceHandle, PropertyType>(*this, name);
    }

    /// This method creates a property about the vertices of the mesh.
    template<typename PropertyType>
    VertexProperty<PropertyType> CreateVertexProperty(const char *name) {
      return OpenMesh::getOrMakeProperty<VertexHandle, PropertyType>(*this, name);
    }

    /// This method creates a property about the edges of the mesh.
    template<typename PropertyType>
    EdgeProperty<PropertyType> CreateEdgeProperty(const char *name) {
      return OpenMesh::getOrMakeProperty<EdgeHandle, PropertyType>(*this, name);
    }

    /// This method creates a temporary property about the faces of the mesh.
    template<typename PropertyType>
    FaceProperty<PropertyType> CreateTemporaryFaceProperty() {
      return OpenMesh::makeTemporaryProperty<FaceHandle, PropertyType>(*this);
    }

    /// This method creates a temporary property about the vertices of the mesh.
    template<typename PropertyType>
    VertexProperty<PropertyType> CreateTemporaryVertexProperty() {
      return OpenMesh::makeTemporaryProperty<VertexHandle, PropertyType>(*this);
    }

    /// This method creates a temporary property about the edges of the mesh.
    template<typename PropertyType>
    EdgeProperty<PropertyType> CreateTemporaryEdgeProperty() {
      return OpenMesh::makeTemporaryProperty<EdgeHandle, PropertyType>(*this);
    }

    /// This method applies a symmetry by a plane of equation z = h.
    void SymmetryHorizontalPlane(const double &height);

    /// This method translates the mesh.
    void Translate(const Vector3d &t);

    /// This method rotates the mesh, based on Cardan angles
    void Rotate(double phi, double theta, double psi);

    //TODO:
//    void AffineTransform(const Eigen::Transform<double, 3, Eigen::Affine>& transform);

    // TODO
//    bool IsWatertight() const;

    /// This method flips the face normals.
    void FlipFaceNormals();

    /// This method flips the vertex normals.
    void FlipVertexNormals();

    /// This method flips the halfegde normals.
    void FlipHalfedgeNormals();

    /// This method flips the normals (face, vertice, halfedge).
    void FlipAllNormals();

    /// Computes triangular faces surface integration of some polynomial integrands using analytical formulas
    /// established by transforming surface integrals into contour integrals and deriving analytical expressions.
    /// Extended from Eberly... https://d-ice.gitlab.host/common/technical_reports/mesh-integrals
    void CalcFacePolynomialIntegrals(const Mesh::FaceHandle &fh);

    /// This method returns the tables of vertices and faces of a mesh.
    void Vertices_Faces(std::vector<Vector3d> &vertices, std::vector<Eigen::VectorXi> &faces);

  };

  /// This function merges two meshes into a single mesh.
  meshoui::Mesh Merge(meshoui::Mesh &mesh_1, meshoui::Mesh &mesh_2);

  /// This function merges several meshes into a single mesh.
  meshoui::Mesh Merge(std::vector<meshoui::Mesh*> meshes);

  meshoui::Mesh operator+(meshoui::Mesh &mesh_1, meshoui::Mesh &mesh_2);

  void Write_OBJ(Mesh &mesh, const std::string &obj_filename);

} // end namespace meshoui

#endif // MESHOUI_LOADER_H

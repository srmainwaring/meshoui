//
// Created by frongere on 18/03/2020.
//

#ifndef MESHOUI_CLIPPER_H
#define MESHOUI_CLIPPER_H

#include <memory>

#include "maths.h"
#include "mesh.h"
#include "plane.h"


namespace meshoui {

  /// Enum flags to locate vertices with respect to clipping surface
  enum VertexPositionWRTClippingSurface {
    // On pourrait du coup plutot utiliser les fonctions d'ajout dynamique de proprietes !!
    VP_ABOVE_SURFACE = 0,
    VP_ON_SURFACE = 1,
    VP_UNDER_SURFACE = 2,
    VP_UNDEFINED = -1
  };

  /// Enum flags to specify how is located a triangle with respect to the clipping surface
  /// Position code is composed of 3 digit AOU where A is the number of vertices above the clipping surface,
  /// O the number of vertices on and U the number of vertices under.
  enum FacePositionType {

    FPT_003 = 3,  // TODO : simplifier la representation de cas en enum, on peut certainemet nommer de maniere unique les cas entierement mouille, sec ou a couper
    //  -----------  totally wet
    //       *
    //      / \
    //     /   \
    //    *-----*
    FPT_012 = 12,
    //   ------*------ totally wet
    //        / \
    //       /   \
    //      *-----*
    FPT_021 = 21,
    //  ----*-----*---- totally wet
    //       \   /
    //        \ /
    //         *
    FPT_030 = 30,
    //  ----*----*----*  Lying on the clipping surface, should be removed

    FPT_102 = 102,
    //         *    Face to clip
    //        / \
    //    ---o---o---
    //      /     \
    //     *-------*
    FPT_111 = 111,
    //          *               *    Face to clip
    //         /|               |\
    //        / |               | \
    //    ---*--o---    or   ---o--*---
    //        \ |               | /
    //         \|               |/
    //          *               *
    FPT_120 = 120,
    //          *  totally dry
    //         / \
    //        /   \
    //   ----*-----*----

    FPT_201 = 201,
    //       *-------*  Face to clip
    //        \     /
    //      ---o---o---
    //          \ /
    //           *
    FPT_210 = 210,
    //        *-----*  totally dry
    //         \   /
    //          \ /
    //       ----*----
    FPT_300 = 300,
    //           *  totally dry
    //          / \
    //         /   \
    //        *-----*
    //
    //   ----------------

    FPT_UNDEFINED = -1
  };

  /// The clipper class which is mostly a functor
  template<class ClippingSurfaceType>
  class Clipper {

   public:

    explicit Clipper(std::shared_ptr<ClippingSurfaceType> clipping_surface) :
        m_clippingSurface(clipping_surface),
        m_Threshold(1e-4),
        m_ProjectionThresholdRatio(1e-3) {}

    /// Returns a clipped copy of the mesh by the clipping surface (does not modifies the input mesh)
    Mesh ClipCopy(const Mesh &mesh);

    /// Clips the input mesh inplace with respect to the clipping surface
    void ClipIt(Mesh &mesh);

    // TODO: remettre en place les capacites de projection qui ont ete desactivees en imposant des valeurs les
    // invalidant dans le constructeur...

//    /// Set the threshold used for crossing and classifying computations
//    /// \param eps threshold
//    void SetThreshold(double eps);
//
//    /// Set the threshold used for projection computations
//    /// \param projectionThresholdRatio threshold
//    void SetProjectionThresholdRatio(double projectionThresholdRatio);

   private:
    /// Initialize the mesh clipper
    void Initialize(Mesh &clipped_mesh);

    /// Clear the mesh
    void Clear(Mesh &clipped_mesh);

    /// This function classify the vertices wrt the clipping surface.
    void ClassifyVertices(Mesh &clipped_mesh);

    /// This function computes the distance wrt the clipping surface and classifies the nodes.
    /// \param vh vertex to be classified
    /// \return vertex position with respect to the clipping surface
    VertexPositionWRTClippingSurface ClassifyVertex(Mesh &clipped_mesh, const Mesh::VertexHandle &vh) const;

    /// This function classfies faces wrt the incident clipping surface.
    /// \param fh face to be classified
    /// \return face position with respect to the clipping surface
    FacePositionType ClassifyFace(Mesh &clipped_mesh, const Mesh::FaceHandle &fh);

    /// Clip the mesh with the given clipping surface
    void Clip(Mesh &clipped_mesh);

    void UpdateModifiedFaceProperties(Mesh &clipped_mesh, Mesh::FaceHandle fh);

    void ClipFace(Mesh &clipped_mesh, const Mesh::FaceHandle &fh);

    bool HasProjection(Mesh &clipped_mesh, const Mesh::FaceHandle &fh);

    /// Performs the clipping of a halfedge, the computation of the intersection node,
    /// the creation of new panels and the deletion of useless panels and vertices.
    void ProcessHalfEdge(Mesh &clipped_mesh, Mesh::HalfedgeHandle heh);

    /// Prepare face for deletion at the end of clipping procedure
    void FlagFaceToBeDeleted(const Mesh::FaceHandle &fh);

    void FlagVertexAdjacentFacesToBeDeleted(Mesh &clipped_mesh, const Mesh::VertexHandle &vh);

    /// Tells if the input edge is crossing the clipping surface
    bool IsEdgeCrossing(Mesh &clipped_mesh, const Mesh::EdgeHandle &eh);

    /// Tells if the input half-edge is crossing the clipping surface
    bool IsHalfEdgeCrossing(Mesh &clipped_mesh, const Mesh::HalfedgeHandle &heh);

    /// Tells if the input half-edge is crossing the clipping surface from above to below
    bool IsHalfEdgeDownCrossing(Mesh &clipped_mesh, const Mesh::HalfedgeHandle &heh);

    /// Tells if the input half-edge is crossing the clipping surface from below to above
    bool IsHalfEdgeUpCrossing(Mesh &clipped_mesh, const Mesh::HalfedgeHandle &heh);

    /// Finds the first half-edge that is crossing the clipping surface from below to above
    Mesh::HalfedgeHandle FindUpcrossingHalfEdge(Mesh &clipped_mesh, const Mesh::FaceHandle &fh);

    /// Finds the first half-edge that is crossing the clipping surface from above to below
    Mesh::HalfedgeHandle FindDowncrossingHalfEdge(Mesh &clipped_mesh, const Mesh::FaceHandle &fh);

    /// Insert a new vertex in the mesh that is intersection between the clipping surface and the input half-edge
    Mesh::VertexHandle InsertIntersectionVertex(Mesh &clipped_mesh, const Mesh::HalfedgeHandle &heh);

//    double GetVertexDistanceToSurface(const Mesh::VertexHandle &vh) const;

    void ApplyFaceDeletion(Mesh &clipped_mesh);

    /// Final cleaning of the mesh (deletion of faces, update of every property, etc.).
    void Finalize(Mesh &clipped_mesh);

   private:
    void AddVertexPositionWRTClippingSurfaceProperty(Mesh &clipped_mesh);


   private:

//    / Initial mesh.
//    Mesh m_mesh;

    /// Clipping surface, by default the plane z = 0.
    std::shared_ptr<ClippingSurfaceType> m_clippingSurface;

    double m_Threshold;
    double m_ProjectionThresholdRatio;

    /// Vector to store the faces which are on and/or above the incident free surface and have to be deleted.
    std::vector<Mesh::FaceHandle> c_FacesToDelete;

    /// Vector to store the faces which need to be clipped.
    std::vector<Mesh::FaceHandle *> c_FacesToUpdate;

  };

}  // end namespace meshoui

// Including the implementation file which is inlined as we are using a templated class
#include "meshoui/clipper.hxx"

#endif //MESHOUI_CLIPPER_H

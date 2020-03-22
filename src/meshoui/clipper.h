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


  /**
  * \class FrClippingSurface
  * \brief Class for dealing with the clipping incident wave field.
  */
  class ClippingSurface {

//   protected:

//    double m_ThresholdDichotomy = 1e-4;     ///< threshold for the dichotomy in the intersection computation

//    Point m_bodyPosition = {0., 0., 0.};   ///< horizontal position of the body, related to the mesh to be clipped

   public:

//    /// Set the body position in world reference frame, for correction in ClippingWaveSurface::GetDistance
//    /// \param bodyPos
//    void SetBodyPosition(Vector3d bodyPos);

    /// This function gives the distance to the clipping surface
    virtual double GetDistance(const Mesh::Point &point) const = 0;

    /// This function gives the intersection node position between an edge and an incident wave field.
    virtual Mesh::Point GetIntersection(const Mesh::Point &p0, const Mesh::Point &p1) = 0;

  };


  /**
  * \class FrClippingPlane
  * \brief Class used when the clipping incident wave field is a HORIZONTAL plane.
  */
  class ClippingPlane : public ClippingSurface {

   private:

    std::unique_ptr<Plane> m_plane;     ///< plane used for clipping

   public:

    ClippingPlane() : m_plane(std::make_unique<Plane>()) {}

//    explicit ClippingPlane(const std::shared_ptr<Plane> &plane);;
//
    /// This function gives the distance to the plane.
    double GetDistance(const Mesh::Point &point) const override;

    /// This function gives the intersection node position between an edge and the plane.
    Mesh::Point GetIntersection(const Mesh::Point &p0, const Mesh::Point &p1) override;
//
//    geom::Plane *GetPlane() const;

  };



//  /**
//  * \class FrClippingWavesSurface
//  * \brief Class for clipping a mesh by an arbitrary incident wave field.
//  */
//  class FrClippingWaveSurface : public FrClippingSurface {
//
//   private:
//
//    FrFreeSurface *m_freeSurface;   ///< free surface used for clipping
//
//   public:
//
//    explicit FrClippingWaveSurface(FrFreeSurface *freeSurface) : m_freeSurface(freeSurface) {};
//
//    /// This function gives the distance to the incident wave.
//    double GetDistance(const FrMesh::Point &point) const override;
//
//    /// This function performs a bisection method to track the intersection node.
//    FrMesh::Point GetIntersection(const FrMesh::Point &p0, const FrMesh::Point &p1) override;
//  };

  enum VertexPositionWRTClippingSurface {
    // On pourrait du coup plutot utiliser les fonctions d'ajout dynamique de proprietes !!
    VP_ABOVE_SURFACE = 0,
    VP_ON_SURFACE = 1,
    VP_UNDER_SURFACE = 2,
    VP_UNDEFINED = -1
  };

  enum FacePositionType {
    // Position code is composed of 3 digit AOU where A is the number of vertices above the clipping surface,
    // O the number of vertices on and U the number of vertices under.
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

  template<typename ClippingSurfaceType>
  class Clipper {

   public:

    Clipper() : m_clippingSurface(std::make_unique<ClippingSurfaceType>()), m_mesh(nullptr) {}

    /// This function gives the clipping surface.
//    ClippingSurface *GetClippingSurface();

    /// Performs the clipping on the specified mesh
    void Apply(Mesh *mesh); // TODO:

    /// Set the threshold used for crossing and classifying computations
    /// \param eps threshold
    void SetThreshold(double eps);

    /// Set the threshold used for projection computations
    /// \param projectionThresholdRatio threshold
    void SetProjectionThresholdRatio(double projectionThresholdRatio);

    /// Set the clipping surface to be used
    /// \param clippingSurface clipping surface
//    void SetClippingSurface(std::shared_ptr<ClippingSurface> clippingSurface);

   private:

    /// Initialize the mesh clipper
    void Initialize();

    /// Clear the mesh
    void Clear();

    /// This function classify the vertices wrt the clipping surface.
    void ClassifyVertices();

    /// This function computes the distance wrt the clipping surface and classifies the nodes.
    /// \param vh vertex to be classified
    /// \return vertex position with respect to the clipping surface
    VertexPositionWRTClippingSurface ClassifyVertex(const Mesh::VertexHandle &vh) const;

    /// This function classfies faces wrt the incident clipping surface.
    /// \param fh face to be classified
    /// \return face position with respect to the clipping surface
    FacePositionType ClassifyFace(const Mesh::FaceHandle &fh);

    /// Clip the mesh with the given clipping surface
    void Clip();

    void UpdateModifiedFaceProperties(FaceHandle fh);

    void ProcessFace(const Mesh::FaceHandle &fh);

    bool HasProjection(const Mesh::FaceHandle &fh);

    void ProcessHalfEdge(Mesh::HalfedgeHandle heh);

    void FlagFaceToBeDeleted(const Mesh::FaceHandle &fh);

    void FlagVertexAdjacentFacesToBeDeleted(const Mesh::VertexHandle &vh);

    bool IsEdgeCrossing(const Mesh::EdgeHandle &eh);

    bool IsHalfEdgeCrossing(const Mesh::HalfedgeHandle &heh);

    bool IsHalfEdgeDownCrossing(const Mesh::HalfedgeHandle &heh);

    bool IsHalfEdgeUpCrossing(const Mesh::HalfedgeHandle &heh);

    Mesh::HalfedgeHandle FindUpcrossingHalfEdge(const Mesh::FaceHandle &fh);

    Mesh::HalfedgeHandle FindDowncrossingHalfEdge(const Mesh::FaceHandle &fh);

    Mesh::VertexHandle InsertIntersectionVertex(const Mesh::HalfedgeHandle &heh);

    double GetVertexDistanceToSurface(const Mesh::VertexHandle &vh) const;

    void ApplyFaceDeletion();

    void Finalize();


   private:

    /// Initial mesh.
    Mesh *m_mesh;

    /// Clipping surface, by default the plane z = 0.
    std::unique_ptr<ClippingSurfaceType> m_clippingSurface;

    double m_Threshold = 1e-4;
    double m_ProjectionThresholdRatio = 1 / 4.;

    /// Vector to store the faces which are on and/or above the incident free surface and have to be deleted.
    std::vector<Mesh::FaceHandle> c_FacesToDelete;

    /// Vector to store the faces which need to be clipped.
    std::vector<Mesh::FaceHandle *> c_FacesToUpdate;

  };


}  // end namespace meshoui

#include "meshoui/clipper.hxx"

#endif //MESHOUI_CLIPPER_H

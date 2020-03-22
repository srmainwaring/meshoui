//
// Created by frongere on 18/03/2020.
//

#include "clipping_surfaces.h"

namespace meshoui {

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::Apply(Mesh *mesh) {

    m_mesh = mesh;

    // Adding dynamic property for vertex positions wrt clipping surface
    AddVertexPositionWRTClippingSurfaceProperty();

    // Partition of the mesh.
    Initialize();

    // Clipping.
    Clip();

    // Clipped mesh cleaning (deletion of faces, update of every property, etc.).
    Finalize();
  }

//  template<typename ClippingSurfaceType>
//  void Clipper<ClippingSurfaceType>::SetProjectionThresholdRatio(double projectionThresholdRatio) {
//    m_ProjectionThresholdRatio = projectionThresholdRatio;
//  }
//
//  template<typename ClippingSurfaceType>
//  void Clipper<ClippingSurfaceType>::SetThreshold(double eps) {
//    m_Threshold = eps;
//  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::Initialize() {

    // Vector to store the faces to delete.
    c_FacesToDelete.reserve(
        m_mesh->n_faces()); // Worst case (fully under the clipping surface) PYW: above the clipping surface?

    // Vector to store the faces to update.
    c_FacesToUpdate.clear();

    // Partition of the vertices of the mesh.
    ClassifyVertices();
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::ClassifyVertices() {

    // Adding a dynamic property to the mesh
    auto vprop = m_mesh->GetVertexProperty<VertexPositionWRTClippingSurface>("vertex_position_wrt_clipping_surface");

    for (auto vh : m_mesh->vertices()) {

      // Computation of the distance to the plane and classification of the nodes.
      auto vPos = ClassifyVertex(vh);

      // Storage of the partition.
      m_mesh->property(*vprop, vh) = vPos; // TODO: verifier

    }
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::Clear() {
    m_mesh->clear();
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::Clip() {
    // Loop over faces.
    for (Mesh::FaceIter fh_iter = m_mesh->faces_begin(); fh_iter != m_mesh->faces_end(); ++fh_iter) {
      ClipFace(*fh_iter);
    }
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::UpdateModifiedFaceProperties(Mesh::FaceHandle fh) {
    Mesh::FFIter ff_iter = m_mesh->ff_iter(fh);
    for (; ff_iter.is_valid(); ++ff_iter) {
      m_mesh->update_normal(*ff_iter); // TODO: remettre en place!  // FIXME: faire un UpdateProperties
      m_mesh->CalcFacePolynomialIntegrals(*ff_iter); //FIXME
    }
  }

  template<typename ClippingSurfaceType>
  bool Clipper<ClippingSurfaceType>::HasProjection(const Mesh::FaceHandle &fh) {

    bool anyVertexProjected = false; // Initialization.
    double dist, edge_length;
    double distMin = 1e10;

    Mesh::HalfedgeHandle oheh;
    Mesh::VertexHandle vh; // Vertex.
    Mesh::Point P0, P1, Pi, Pi_final; // Vertices.
    Mesh::VertexOHalfedgeIter voh_iter; // Vertex outgoing halfedge iterator.

    // Iterating on vertices of the face to clip.
    Mesh::FaceVertexIter fv_iter = m_mesh->fv_iter(fh);

    for (; fv_iter.is_valid(); ++fv_iter) {

      // First vertex of the edge.
      vh = *fv_iter;
      P0 = m_mesh->point(vh);

      // Iterating on outgoing halfedges to get the shortest edge path.
      // A vertex has several outgoing halfedges, pointing to all its neighbooring vertices.
      voh_iter = m_mesh->voh_iter(vh);
      for (; voh_iter.is_valid(); ++voh_iter) {

        // outgoing halfedge
        oheh = *voh_iter;

        // Is the halfedge clipping the surface?
        if (IsHalfEdgeCrossing(oheh)) {

          // Second vertex of the edge.
          P1 = m_mesh->point(m_mesh->to_vertex_handle(oheh));

          // Get the intersection node position.

          Pi = m_clippingSurface->GetIntersection(P0, P1);

          dist = (Pi - P0).norm();
          edge_length = (P1 - P0).norm(); // TODO: utiliser le precalcul...

          // If the intersection node is too close from the vertex vh, it will be projected.
          if (dist < edge_length * m_ProjectionThresholdRatio) {
            if (dist < distMin) {

              // Storage of the intersection node which is the closest of the vertex vh.
              Pi_final = Pi;
              distMin = dist;
            }
            anyVertexProjected = true;
          }
        }
      }

      auto vprop = m_mesh->GetVertexProperty<VertexPositionWRTClippingSurface>("vertex_position_wrt_clipping_surface");

      // If anyVertexProjected is true, the vertex vh will be projected to one of its neighboring vertex.
      // The choice is made based on the minimum distance to them, by using distMin.
      // If the vertex is projected, no new vertex will be created.
      if (anyVertexProjected) {
        // Placing vh on Pi.
        m_mesh->point(vh) = Pi_final;
        m_mesh->property(*vprop, vh) = VP_ON_SURFACE; // TODO: verifier

        break; // We anyVertexProjected only one vertex per face at a time as other will be processed by adjacent faces
      }
    }

    return anyVertexProjected;
  }

  template<typename ClippingSurfaceType>
  VertexPositionWRTClippingSurface Clipper<ClippingSurfaceType>::ClassifyVertex(const Mesh::VertexHandle &vh) const {
    double distance = m_clippingSurface->GetDistance(m_mesh->point(vh));
    // TODO: On peut projeter le vertex si distance est petit (si plus petit que meanEdgeLength * m_threshold)

    VertexPositionWRTClippingSurface vPos;

    if (fabs(distance) < m_Threshold) {
      vPos = VP_ON_SURFACE;
    } else if (distance > 0.) {
      vPos = VP_ABOVE_SURFACE;
    } else {
      vPos = VP_UNDER_SURFACE;
    }

    return vPos;
  }

  template<typename ClippingSurfaceType>
  FacePositionType
  Clipper<ClippingSurfaceType>::ClassifyFace(const Mesh::FaceHandle &fh) { // TODO: renvoyer le resultat !!

    FacePositionType fPos;
    unsigned int nbAbove, nbUnder;
    nbAbove = nbUnder = 0;

    auto vprop = m_mesh->GetVertexProperty<VertexPositionWRTClippingSurface>("vertex_position_wrt_clipping_surface");

    // Counting the number of vertices above and under the clipping surface.
    Mesh::FaceVertexIter fv_iter = m_mesh->fv_iter(fh);
    for (; fv_iter.is_valid(); ++fv_iter) {
      auto vpos = m_mesh->property(*vprop, *fv_iter);

      if (vpos == VP_ABOVE_SURFACE) {
        ++nbAbove;
      }

      if (vpos == VP_UNDER_SURFACE) {
        ++nbUnder;
      }
    }

    // Face under or on the clipping surface.
    if (nbAbove == 0) {
      if (nbUnder == 0) {
        fPos = FPT_030;
      } else if (nbUnder == 1) {
        fPos = FPT_021;
      } else if (nbUnder == 2) {
        fPos = FPT_012;
      } else {
        fPos = FPT_003;
      }
      // Face crossing the clipping surface.
    } else if (nbAbove == 1) {
      if (nbUnder == 0) {
        fPos = FPT_120;
      } else if (nbUnder == 1) {
        fPos = FPT_111;
      } else {
        fPos = FPT_102;
      }
      // Face crossing the clipping surface.
    } else if (nbAbove == 2) {
      if (nbUnder == 0) {
        fPos = FPT_210;
      } else {
        fPos = FPT_201;
      }
      // Face above the clipping surface.
    } else {
      fPos = FPT_300;
    }

    return fPos;
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::ClipFace(const Mesh::FaceHandle &fh) {

    Mesh::HalfedgeHandle heh_down, heh_up;
    Mesh::FaceEdgeIter fe_iter;

    // Classification of faces wrt the incident wave field.
    FacePositionType fPos = ClassifyFace(fh);

    switch (fPos) { // Pourquoi ne pas faire tout ici ?
      case FPT_003:
        // Totally wet, we keep the face
        break;
      case FPT_012:
        // Totally wet, we keep the face
        break;
      case FPT_021:
        // Totally wet, we keep the face
        break;
      case FPT_030:
        // Face on the clipping surface.
        FlagFaceToBeDeleted(fh);
        break;
      case FPT_102:
//                    c_FacesToUpdate.emplace_back(new FaceHandle(const_cast<FaceHandle&>(fh)));

        // If HasProjection = True, some vertices are projected to the intersection nodes and no new face or vertex is added.

        if (HasProjection(fh)) {
          // The function ClipFace is run again to update the classfication?
          ClipFace(fh);
        } else {

          // Clipping the panel, creation of new nodes on the intersection curve and new faces.
          // TODO: avoir une methode qui a la fois ajouter le vertex et plitte l'edge...

          // Unique half edge oriented downward the intersection curve.
          heh_down = FindDowncrossingHalfEdge(fh);

          // Unique half edge oriented upward the intersection curve.
          heh_up = FindUpcrossingHalfEdge(fh);

          // Clipping of the upward and downward harfedges, creation of new vertices and faces and deletion of the useless ones.
          ProcessHalfEdge(heh_down);
          ProcessHalfEdge(heh_up);

        }
//                    UpdateModifiedFaceProperties(fh);
        break;

      case FPT_111:
//                    c_FacesToUpdate.emplace_back(new FaceHandle(const_cast<FaceHandle&>(fh)));

        if (HasProjection(fh)) {
          // The function ProcessFase is run again to update the classfication?
          ClipFace(fh);
        } else {
          fe_iter = m_mesh->fe_iter(fh);
          for (; fe_iter.is_valid(); ++fe_iter) {
            if (IsEdgeCrossing(*fe_iter)) {
              break;
            }
          }

          ProcessHalfEdge(m_mesh->halfedge_handle(*fe_iter, 0));
        }
//                    UpdateModifiedFaceProperties(fh);
        break;

      case FPT_120:
        // Face above and on the clipping surface.
        FlagFaceToBeDeleted(fh);
        break;

      case FPT_201:
//                    c_FacesToUpdate.emplace_back(new FaceHandle(const_cast<FaceHandle&>(fh)));

        if (HasProjection(fh)) {
          // The function ProcessFase is run again to update the classfication?
          ClipFace(fh);
        } else {

          // TODO: avoir une methode qui a la fois ajouter le vertex et plitte l'edge...
          heh_down = FindDowncrossingHalfEdge(fh);
          heh_up = FindUpcrossingHalfEdge(fh);

          ProcessHalfEdge(heh_down);
          ProcessHalfEdge(heh_up);
        }
//                    UpdateModifiedFaceProperties(fh);
        break;

      case FPT_210:
        // Face above and on the clipping surface.
        FlagFaceToBeDeleted(fh);
        break;

      case FPT_300:
        // Face above the clipping surface.
        FlagFaceToBeDeleted(fh);
        break;

      case FPT_UNDEFINED:
        // TODO: throw exception
        break;

    }
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::ProcessHalfEdge(Mesh::HalfedgeHandle heh) {

    Mesh::VertexHandle vh;  // FIXME: au final, on va juste effectuer l'intersection vu qu'on va projeter les

    // Intersection node.
    vh = InsertIntersectionVertex(heh);

    // Clipping, updating of the mesh, deletion of useless panels and vertices.
    m_mesh->split(m_mesh->edge_handle(heh), vh);

    // Updating the faces to delete.
    FlagVertexAdjacentFacesToBeDeleted(vh);

  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::FlagFaceToBeDeleted(const Mesh::FaceHandle &fh) {
    c_FacesToDelete.push_back(fh);
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::FlagVertexAdjacentFacesToBeDeleted(const Mesh::VertexHandle &vh) {

    Mesh::VertexFaceIter vf_iter = m_mesh->vf_iter(vh);
    FacePositionType fPos;
    for (; vf_iter.is_valid(); ++vf_iter) {
      fPos = ClassifyFace(*vf_iter);
      if (fPos == FPT_120 || fPos == FPT_210 || fPos == FPT_030) {
        FlagFaceToBeDeleted(*vf_iter);
      }
    }
  }

  template<typename ClippingSurfaceType>
  bool Clipper<ClippingSurfaceType>::IsEdgeCrossing(const Mesh::EdgeHandle &eh) {
    Mesh::HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);

    double dz_0 = m_clippingSurface->GetDistance(m_mesh->point(m_mesh->from_vertex_handle(heh)));
    double dz_1 = m_clippingSurface->GetDistance(m_mesh->point(m_mesh->to_vertex_handle(heh)));
    double prod = dz_0 * dz_1;
    bool out;

    if (fabs(dz_0) < m_Threshold || fabs(dz_1) < m_Threshold) {
      out = false;
    } else {
      // If prod is negative, the two vertices are not on the same side of the clipping surface.
      out = (prod < 0.);
    }

    return out;
  }

  template<typename ClippingSurfaceType>
  bool Clipper<ClippingSurfaceType>::IsHalfEdgeCrossing(const Mesh::HalfedgeHandle &heh) {
    return IsEdgeCrossing(m_mesh->edge_handle(heh));
  }

  template<typename ClippingSurfaceType>
  bool Clipper<ClippingSurfaceType>::IsHalfEdgeDownCrossing(const Mesh::HalfedgeHandle &heh) {

    double dz_from = m_clippingSurface->GetDistance(m_mesh->point(m_mesh->from_vertex_handle(heh)));
    double dz_to = m_clippingSurface->GetDistance(m_mesh->point(m_mesh->to_vertex_handle(heh)));
    bool out;
    if (fabs(dz_from) < m_Threshold || fabs(dz_to) < m_Threshold) {
      out = false;
    } else {
      out = (dz_from > 0. && dz_to < 0.);
    }
    return out;
  }

  template<typename ClippingSurfaceType>
  bool Clipper<ClippingSurfaceType>::IsHalfEdgeUpCrossing(const Mesh::HalfedgeHandle &heh) {
    // FIXME: cette methode est boilerplate par rapport a sa version down crossing !!!

    double dz_from = m_clippingSurface->GetDistance(m_mesh->point(m_mesh->from_vertex_handle(heh)));
    double dz_to = m_clippingSurface->GetDistance(m_mesh->point(m_mesh->to_vertex_handle(heh)));
    bool out;
    if (fabs(dz_from) < m_Threshold || fabs(dz_to) < m_Threshold) {
      out = false;
    } else {
      out = (dz_from < 0. && dz_to > 0.);
    }
    return out;
  }

  template<typename ClippingSurfaceType>
  Mesh::HalfedgeHandle Clipper<ClippingSurfaceType>::FindUpcrossingHalfEdge(const Mesh::FaceHandle &fh) {

    // TODO: throw error if no upcrossing halfedge is found
    Mesh::HalfedgeHandle heh = m_mesh->halfedge_handle(fh);
    while (!IsHalfEdgeUpCrossing(heh)) {
      heh = m_mesh->next_halfedge_handle(heh);
    }
    return heh;
  }

  template<typename ClippingSurfaceType>
  Mesh::HalfedgeHandle Clipper<ClippingSurfaceType>::FindDowncrossingHalfEdge(const Mesh::FaceHandle &fh) {
    // TODO: throw error if no downcrossing halfedge is found
    Mesh::HalfedgeHandle heh = m_mesh->halfedge_handle(fh);
    unsigned int i = 0;
    while (!IsHalfEdgeDownCrossing(heh) && i < 2) { // TODO: abandonner le i pour le garde fou... --> erreur
      i++;
      heh = m_mesh->next_halfedge_handle(heh);
    }
    return heh;
  }

  template<typename ClippingSurfaceType>
  Mesh::VertexHandle Clipper<ClippingSurfaceType>::InsertIntersectionVertex(const Mesh::HalfedgeHandle &heh) {
    Mesh::Point p0 = m_mesh->point(m_mesh->from_vertex_handle(heh));
    Mesh::Point p1 = m_mesh->point(m_mesh->to_vertex_handle(heh));

    Mesh::Point p_intersection = m_clippingSurface->GetIntersection(
        m_mesh->point(m_mesh->from_vertex_handle(heh)),
        m_mesh->point(m_mesh->to_vertex_handle(heh))
    );

    Mesh::VertexHandle vh = m_mesh->add_vertex(p_intersection);

    auto vprop = m_mesh->GetVertexProperty<VertexPositionWRTClippingSurface>("vertex_position_wrt_clipping_surface");

//    m_mesh->data(vh).SetOn(); // Vertex has been built on the clipping surface
    m_mesh->property(*vprop, vh) = VP_ON_SURFACE;

    return vh;
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::ApplyFaceDeletion() {
    for (Mesh::FaceHandle fh : c_FacesToDelete) {
      if (!m_mesh->status(fh).deleted()) {
        m_mesh->delete_face(fh);
      }
    }
    c_FacesToDelete.clear();
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::Finalize() {
    ApplyFaceDeletion();

    // FIXME: ne fonctionne pas, il faut que ces vecteurs soient initialises avec des elements a tracker
    std::vector<Mesh::VertexHandle *> vh_to_update;
    std::vector<Mesh::HalfedgeHandle *> hh_to_update;
//            std::vector<FaceHandle*> fh_to_update;

    m_mesh->garbage_collection(vh_to_update, hh_to_update, c_FacesToUpdate);

    m_mesh->UpdateAllProperties(); // FIXME: il ne faudrait mettre a jour que les pptes de facettes ayant ete mofifiees ou nouvellement crees

    // We have to update the modified and new faces properties
//            DMesh::FaceFaceIter ff_iter;
//            for (FaceHandle* fh_ptr : c_FacesToUpdate) {
//                if (!fh_ptr->is_valid()) continue;
//
//                m_mesh->update_normal(*fh_ptr);
//                // TODO: faire l'update aussi des centres !!
//                m_mesh->CalcFacePolynomialIntegrals(*fh_ptr);
//
//                // Updating the surrounding faces too
//                ff_iter = m_mesh->ff_iter(*fh_ptr);
//                for (; ff_iter.is_valid(); ++ff_iter) {
//                    m_mesh->update_normal(*ff_iter);
//                    m_mesh->CalcFacePolynomialIntegrals(*ff_iter);
//                }
//
//
//                // Updating also the surrounding faces
//            }
//
//            // Updating the modified and new vertices properties
//            for (VertexHandle* vh_ptr : vh_to_update) {
//                m_mesh->update_normal(*vh_ptr);
//            }

    // FIXME: il faut egalement detruire les pptes dynamiques ajoutees temporairement au maillage
//    m_mesh->property


    m_mesh = nullptr;
  }

  template<typename ClippingSurfaceType>
  void Clipper<ClippingSurfaceType>::AddVertexPositionWRTClippingSurfaceProperty() {
    // Dynamic property creation
    auto vprop = m_mesh->CreateVertexProperty<VertexPositionWRTClippingSurface>("vertex_position_wrt_clipping_surface");

    // Setting every vertex to undefined state for this property
    for (const auto &vh : m_mesh->vertices()) {
      m_mesh->property(*vprop, vh) = VP_UNDEFINED;
    }
  }

//  template<typename ClippingSurfaceType>
//  Planar3DPolygon Clipper<ClippingSurfaceType>::ExtractClippedPolygonSet() const {
//
//    auto clipping_plane = dynamic_cast<ClippingPlane*>(m_clippingSurface.get());
//
//    if (!clipping_plane) {
//      std::cerr << "Cannot extract a planar polygon set if the clipping surface is not a plane" << std::endl;
//      exit(EXIT_FAILURE);
//    }
//
////  TODO: faire des tests pour verifier que le polygone frontiere suive bien la surface
//    // on peut avoir une methode qui accepte une Clipping Surface et qui verifie qu'on a bien une intersection
//    // avec la surface de decoupe.
//
//    Planar3DPolygonSet polygon_set(clipping_plane->GetPlane());
//
//    // buffer to keep track of the visited halfedges to reinit the tag after it is being used
//    std::vector<Mesh::HalfedgeHandle> tagged_halfedges;
//
//    Mesh::HalfedgeHandle heh_init, heh;
//
//    heh_init = FindFirstUntaggedBoundaryHalfedge();
////    while (heh_init.idx() != -1) {  // TODO : voir s'il n'y a pas de methode heh_init.is_valid()
////      Polygon polygon;
////      std::vector<Position> vertexList;
////
////      polygon.push_back(heh_init);
////      status(heh_init).set_tagged(true);
////      tagged_halfedges.push_back(heh_init);
////
////      auto vertex = OpenMeshPointToVector3d<Position>(point(from_vertex_handle(heh_init)));
////      vertexList.push_back(vertex);
////      vertex = OpenMeshPointToVector3d<Position>(point(to_vertex_handle(heh_init)));
////      vertexList.push_back(vertex);
////
////      // Circulating over the boundary from heh_init until the polygon is closed
////      heh = next_halfedge_handle(heh_init);
////      while (heh != heh_init) {
////        polygon.push_back(heh);
////        status(heh).set_tagged(true);
////        tagged_halfedges.push_back(heh);
////        vertex = OpenMeshPointToVector3d<Position>(point(to_vertex_handle(heh)));
////        vertexList.push_back(vertex);
////        heh = next_halfedge_handle(heh);
////      }
////
////      polygonSet.push_back(FrPolygon(vertexList, NWU));
////
////      heh_init = FindFirstUntaggedBoundaryHalfedge();
////
////    }
////
////    // Removing tags
////    for (HalfedgeHandle heh_ptr : tagged_halfedges) {
////      status(heh_ptr).set_tagged(false);
////    }
////
////    m_polygonSet = polygonSet;
//
//  }


}  // end namespace meshoui

// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_INTEGRATION_H
#define MESHOUI_INTEGRATION_H

#include "maths.h"
#include "mesh.h"
#include "MathUtils/Integration2dTriangle.h"
#include "IntegrandOnFace.h"

namespace meshoui {
  /**
  * Class for dealing integration over the faces.
  */
  template<typename T>
  class Integration {

   public:

    /// Constructor of the class.
    Integration(IntegrandOnFace<T> *F, const int &order, meshoui::Mesh *mesh)
        : m_integrator(F, order), m_order(order), m_mesh(mesh) {
    }

    /// This function computes the surface integration over a face.
    T Compute(const FaceHandle &fh) {

      // Vertices of the triangle.
      auto heh = m_mesh->halfedge_handle(fh);
      Vector3d P0 = m_mesh->point(m_mesh->from_vertex_handle(heh));

      heh = m_mesh->next_halfedge_handle(heh);
      Vector3d P1 = m_mesh->point(m_mesh->from_vertex_handle(heh));

      heh = m_mesh->next_halfedge_handle(heh);
      Vector3d P2 = m_mesh->point(m_mesh->from_vertex_handle(heh));

      // Computation of the integral.
      T result = m_integrator.Compute(P0, P1, P2);

      return (result);

    }

   private:

    /// Function to be integrated
    mathutils::Integration2dTriangle<T> m_integrator;

    /// Order of the quadrature.
    int m_order;

    /// Meshoui mesh.
    meshoui::Mesh *m_mesh;

  };

}

#endif //MESHOUI_INTEGRATION_H

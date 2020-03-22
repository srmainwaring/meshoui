//
// Created by frongere on 22/03/2020.
//

#include "clipping_surfaces.h"


namespace meshoui {

  double ClippingPlane::GetDistance(const Mesh::Point &point) const {
    return m_plane->GetSignedDistanceToPoint(point);
  }

  Mesh::Point ClippingPlane::GetIntersection(const Mesh::Point &p0, const Mesh::Point &p1) {
    return m_plane->GetIntersectionWithLine(p0, p1);
  }


}  // end namespace meshoui

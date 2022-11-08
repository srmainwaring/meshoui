//
// Created by frongere on 22/03/2020.
//

#ifndef MESHOUI_CLIPPING_SURFACES_H
#define MESHOUI_CLIPPING_SURFACES_H

#include "mesh.h"

namespace meshoui {

/**
  * \class FrClippingSurface
  * \brief Class for dealing with the clipping incident wave field.
  */
  class ClippingSurface {

   public:
    virtual ~ClippingSurface() {}

    /// Gives the distance to the clipping surface
    virtual double GetDistance(const Mesh::Point &point) const = 0;

    /// Gives the intersection point between the segment {p0, p1} and the clipping surface
    virtual Mesh::Point GetIntersection(const Mesh::Point &p0, const Mesh::Point &p1) = 0;

  };


  /**
  * \class FrClippingPlane
  * \brief Class used when the clipping incident wave field is a HORIZONTAL plane.
  */
  class ClippingPlane : public ClippingSurface {

   private:

    std::shared_ptr<Plane> m_plane;     ///< plane used for clipping

   public:
    virtual ~ClippingPlane() {}

    ClippingPlane() : m_plane(std::make_shared<Plane>()) {}

    ClippingPlane(std::shared_ptr<Plane> plane) : m_plane(plane) {}

    // FIXME: les 2 methodes suivantes doivent reposer sur les methodes fournies par les objets geometriques

    /// This function gives the distance to the plane.
    double GetDistance(const Mesh::Point &point) const override;

    /// This function gives the intersection node position between an edge and the plane.
    Mesh::Point GetIntersection(const Mesh::Point &p0, const Mesh::Point &p1) override;

    std::shared_ptr<Plane> GetPlane();

  };

}  // end namespace meshoui



#endif //MESHOUI_CLIPPING_SURFACES_H

//
// Created by frongere on 20/03/2020.
//

#ifndef MESHOUI_PLANE_H
#define MESHOUI_PLANE_H

#include "maths.h"

#include "Eigen/Geometry"

// TODO: ce code devrait se trouver dans mathutils::geometry

namespace meshoui {

  double GetAngleBetweenVectors(const Vector3d &v1, const Vector3d &v2);

  /**
     * Class defining geometrically a plane, based on a point, considered as the plane origin and a normal vector
     */
  class Plane {

   private:
    using HyperPlane = Eigen::Hyperplane<double, 3>;
    using Transform = Eigen::Transform<double, 3, Eigen::Affine>;

   public:
    Plane();

    Plane(const Vector3d &origin, const Vector3d &normal);

    Vector3d origin() const;

    Vector3d normal() const;

    void Set(const Vector3d& origin, const Vector3d& normal);

    void SetOrigin(const Vector3d& origin);

    void SetNormal(const Vector3d& normal);

    double GetSignedDistanceToPoint(const Vector3d &point) const;

    Vector3d NormalProjectPointOnPlane(const Vector3d &point) const;

    double GetPlaneOffset() const;

    bool IsPointOnPlane(const Vector3d &point) const;

    Vector3d GetPointPositionInPlaneLocalFrame(const Vector3d &point);

    /// Point must be on the plane
    Vector2d Get2DPointPositionInPlaneLocalFrame(const Vector3d& point);

    Vector3d GetIntersectionWithLine(const Vector3d& p0, const Vector3d& p1) const;

//   private:
    Vector3d GetPlaneFrameOrigin() const;

    void BuildTransform();

    void Update();

   private:
    HyperPlane m_hyperplane;
    Vector3d m_origin;
    Transform c_transform;
  };

}  // end namespace meshoui



#endif //MESHOUI_PLANE_H

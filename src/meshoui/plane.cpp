//
// Created by frongere on 20/03/2020.
//

#include "plane.h"


namespace meshoui {


  double GetAngleBetweenVectors(const Vector3d &v1, const Vector3d &v2) { // TODO: mettre dans mathutils
    return std::atan2(v1.cross(v2).norm(), v1.dot(v2));
  }

  Plane::Plane() : m_origin(0., 0., 0.), m_hyperplane({0., 0., 1.}, 0.) {
    BuildTransform();
  }

  Plane::Plane(const Vector3d &origin, const Vector3d &normal) :
      m_origin(origin),
      m_hyperplane(normal, origin.dot(normal)) {  // FIXME: retirer le - si pas besoin
    BuildTransform();
  }

  Vector3d Plane::origin() const { return m_origin; }

  Vector3d Plane::normal() const { return m_hyperplane.normal(); }

  double Plane::GetSignedDistanceToPoint(const Vector3d &point) const {
    return m_hyperplane.signedDistance(point) - 2 * m_hyperplane.offset();
    // FIXME: the -2*offset seems to be due to an Eigen BUG into signedDistance method of HyperPlane
  }

  Vector3d Plane::NormalProjectPointOnPlane(const Vector3d &point) const {
    return -m_hyperplane.projection(point);
    // FIXME: see what wrong with Eigen and why we have to negate... Eigen BUG ?
  }

  double Plane::GetPlaneOffset() const {
    return m_hyperplane.offset();
  }

  bool Plane::IsPointOnPlane(const Vector3d &point) const {
    return mathutils::IsClose(GetSignedDistanceToPoint(point), 0.);
  }

  Vector3d Plane::GetPointPositionInPlaneLocalFrame(const Vector3d &point) {
    return c_transform * point;
  }

  Vector2d Plane::Get2DPointPositionInPlaneLocalFrame(const Vector3d &point) {
    if (!IsPointOnPlane(point)) {
      std::cerr << "Point is not on plane" << std::endl;
//        exit(EXIT_FAILURE);
    }
    Vector3d projection = GetPointPositionInPlaneLocalFrame(point);
    assert(mathutils::IsClose(projection.z(), 0.)); // FIXME: redondant avec le if dessus
    return projection.head(2);
  }

  Vector3d Plane::GetIntersectionWithLine(const Vector3d &p0, const Vector3d &p1) const {
    // TODO: utiliser code d'erreur
    auto normal = m_hyperplane.normal();

    Vector3d p0p1 = p1 - p0;
    if (mathutils::IsClose(p0p1.norm(), 0.)) {
      std::cerr << "Points defining the line have the same position" << std::endl;
      return {0., 0., 0.};
    }
    if (mathutils::IsClose(p0p1.dot(normal), 0.)) {
      std::cerr << "Line is parallel to the plane" << std::endl;
      return {0., 0., 0.};
    }

    return p0 + p0p1 * (m_origin - p0).dot(normal) / p0p1.dot(normal);

  }

  Vector3d Plane::GetPlaneFrameOrigin() const {
    return NormalProjectPointOnPlane({0., 0., 0.});
  }

  void Plane::BuildTransform() {
    Vector3d normal = m_hyperplane.normal();
    Vector3d z_axis = {0., 0., 1.};

    Transform transform(Eigen::AngleAxisd(GetAngleBetweenVectors({0., 0., 1.}, normal),
                                          (z_axis.cross(normal)).normalized()));
    transform.translation() = m_origin;
    c_transform = transform.inverse();

  }
}  // end namespace meshoui

//
// Created by frongere on 20/03/2020.
//

#include "polygon.h"


namespace meshoui {

  PolygonSurfaceIntegrals::PolygonSurfaceIntegrals() :
      m_Int_1(0.),
      m_Int_x(0.),
      m_Int_y(0.),
      m_Int_xy(0.),
      m_Int_x2(0.),
      m_Int_y2(0.),
      m_Int_x2y(0.),
      m_Int_y2x(0.),
      m_Int_x3(0.),
      m_Int_y3(0.) {}

  PolygonSurfaceIntegrals::PolygonSurfaceIntegrals(double Int_1, double Int_x, double Int_y, double Int_xy,
                                                   double Int_x2, double Int_y2, double Int_x2y, double Int_y2x,
                                                   double Int_x3, double Int_y3) :
      m_Int_1(Int_1),
      m_Int_x(Int_x),
      m_Int_y(Int_y),
      m_Int_xy(Int_xy),
      m_Int_x2(Int_x2),
      m_Int_y2(Int_y2),
      m_Int_x2y(Int_x2y),
      m_Int_y2x(Int_y2x),
      m_Int_x3(Int_x3),
      m_Int_y3(Int_y3) {}


  double PolygonSurfaceIntegrals::GetSurfaceIntegral(PolygonSurfaceIntegrals::IntegrandType type) const {
    switch (type) {
      case POLY_1:
        return m_Int_1;
      case POLY_X:
        return m_Int_x;
      case POLY_Y:
        return m_Int_y;
      case POLY_XY:
        return m_Int_xy;
      case POLY_X2:
        return m_Int_x2;
      case POLY_Y2:
        return m_Int_y2;
      case POLY_X2Y:
        return m_Int_x2y;
      case POLY_Y2X:
        return m_Int_y2x;
      case POLY_X3:
        return m_Int_x3;
      case POLY_Y3:
        return m_Int_y3;
      default:
        std::cerr << "No integration rule for integrand of type " << type << " for polygons" << std::endl;
        exit(1);
    }
  }

  Planar2DPolygon::Planar2DPolygon() : c_is_up_to_date(false) {}

  void Planar2DPolygon::Add(const Vector2d &point) {
    m_points.push_back(point);
    c_is_up_to_date = false;
  }

  void Planar2DPolygon::Update() const {
    UpdateSurfaceIntegrals();
  }

  void Planar2DPolygon::UpdateSurfaceIntegrals() const {

    double Int1, IntX, IntY, IntXY, IntX2, IntY2, IntX2Y, IntY2X, IntX3, IntY3;
    Int1 = IntX = IntY = IntXY = IntX2 = IntY2 = IntX2Y = IntY2X = IntX3 = IntY3 = 0.;

    Vector2d p0, p1;
    double x0, x1, y0, y1;
    double dx, dy, px, py, a, b;

    p0 = m_points[0];

    for (int i = 1; i < m_points.size(); i++) {

      p1 = m_points[i];
      x0 = p0[0];
      y0 = p0[1];

      x1 = p1[0];
      y1 = p1[1];

      dx = x1 - x0;
      dy = y1 - y0;
      px = x0 + x1;
      py = y0 + y1;
      a = x0 * x0 + x1 * x1;
      b = y0 * y0 + y1 * y1;

      Int1 += dy * px;
      IntX += dy * (px * px - x0 * x1);
      IntY += dx * (py * py - y0 * y1);
//        IntXY += dy * (py * a + 2 * px * (x0 * y0 + x1 * y1)); missing terms from FR analytical dev.
      IntXY += dy * (py * px * px + 2 * (y0 * x0 * x0 + y1 * x1 * x1));
      IntX2 += dy * a * px;
      IntY2 += dx * b * py;
      IntX2Y += dy *
                (py * std::pow(px, 3.) + 3. * std::pow(x0, 3.) * y0 + 3. * std::pow(x1, 3.) * y1 - x0 * x0 * x1 * y1 -
                 x0 * x1 * x1 * y0);
      IntY2X += dx *
                (px * std::pow(py, 3.) + 3. * std::pow(y0, 3.) * x0 + 3. * std::pow(y1, 3.) * x1 - y0 * y0 * y1 * x1 -
                 y0 * y1 * y1 * x0);
      IntX3 += dy * (std::pow(x0, 4.) + std::pow(x0, 3.) * x1 + x0 * std::pow(x1, 3.) + std::pow(x1, 4.));
      IntY3 += dx * (std::pow(y0, 4.) + std::pow(y0, 3.) * y1 + y0 * std::pow(y1, 3.) + std::pow(y1, 4.));

      p0 = p1;
    }

    Int1 /= 2.;
    IntX /= 6.;
    IntY /= -6.;
    IntXY /= 24.;
    IntX2 /= 12.;
    IntY2 /= -12.;
    IntX2Y /= 60.;
    IntY2X /= -60.;
    IntX3 /= 20.;
    IntY3 /= -20.;

    c_surface_integrals = PolygonSurfaceIntegrals(Int1, IntX, IntY, IntXY, IntX2, IntY2, IntX2Y, IntY2X, IntX3, IntY3);
  }

  double Planar2DPolygon::GetSurfaceIntegral(PolygonSurfaceIntegrals::IntegrandType type) const {
    if (!c_is_up_to_date) {
      Update();
      c_is_up_to_date = true;
    }
    return c_surface_integrals.GetSurfaceIntegral(type);
  }


  Planar3DPolygon::Planar3DPolygon(std::shared_ptr<Plane> plane) :
      c_is_up_to_date(false),
      m_plane(plane) {}

  void Planar3DPolygon::AddPoint(const Vector3d &point) {
    if (!m_plane->IsPointOnPlane(point)) {
      std::cerr << "Attempting to add a point which is not in the plane" << std::endl;
      exit(EXIT_FAILURE);
    }
    m_points.push_back(point);
  }

  std::shared_ptr<Plane> Planar3DPolygon::GetPlane() const {
    return m_plane;
  }

  Planar2DPolygon Planar3DPolygon::Get2DPolygon() const {
    if (!c_is_up_to_date) Generate2DPolygon();
    return c_2D_polygon;
  }

  void Planar3DPolygon::Generate2DPolygon() const {
    Planar2DPolygon polygon2D;
    for (const auto &point : m_points) {
      polygon2D.Add(m_plane->Get2DPointPositionInPlaneLocalFrame(point));
    }
    c_2D_polygon = polygon2D;
  }

  void Planar3DPolygon::Update() const {
    Generate2DPolygon();
    c_2D_polygon.Update();
  }

//  AABB2D Planar3DPolygon::GetBoundingBox() const {
//    AABB2D bbox;
//
//    //FIXME : plutôt utiliser GetVertexInPlane? --> semi-OOBB?
//
//    bbox.xmin = m_points[0].x();
//    bbox.xmax = m_points[0].x();
//    bbox.ymin = m_points[0].y();
//    bbox.ymax = m_points[0].y();
//
//    for (auto &vertex : m_points) {
//      bbox.xmin = fmin(bbox.xmin, vertex.x());
//      bbox.xmax = fmax(bbox.xmax, vertex.x());
//      bbox.ymin = fmin(bbox.ymin, vertex.y());
//      bbox.ymax = fmax(bbox.ymax, vertex.y());
//    }
//    return bbox;
//  }



}  // end namespace meshoui

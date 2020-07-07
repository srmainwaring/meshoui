//
// Created by frongere on 20/03/2020.
//

#ifndef MESHOUI_POLYGON_H
#define MESHOUI_POLYGON_H

#include <vector>
#include <memory>

#include "maths.h"
#include "plane.h"


// TODO: placer dans mathutils !!

namespace meshoui {

  class PolygonSurfaceIntegrals {

   public:

    enum IntegrandType {
      // TODO: voir http://www.drdobbs.com/when-enum-just-isnt-enough-enumeration-c/184403955 pour une meilleure
      // gestion des enums
      UNDEFINED_INTEGRAND,
      POLY_1,
      POLY_X,
      POLY_Y,
      POLY_Z,
      POLY_YZ,
      POLY_XZ,
      POLY_XY,
      POLY_X2,
      POLY_Y2,
      POLY_Z2,
      POLY_X3,
      POLY_Y3,
      POLY_Z3,
      POLY_X2Y,
      POLY_Y2Z,
      POLY_Z2X,
      POLY_Y2X,
//            INFINITE_DEPTH_GREEN_FUNCTION, // TODO : en parler avec Camille et Lucas et PYW
//            FINITE_DEPTH_GREEN_FUNCTION
    };

    PolygonSurfaceIntegrals();

    PolygonSurfaceIntegrals(double Int_1, double Int_x, double Int_y, double Int_xy, double Int_x2, double Int_y2,
                            double Int_x2y, double Int_y2x, double Int_x3, double Int_y3);

    /// This function gives the surface integral of a mesh.
    double GetSurfaceIntegral(IntegrandType type) const;

   private:
    double m_Int_1;
    double m_Int_x;
    double m_Int_y;
    double m_Int_xy;
    double m_Int_x2;
    double m_Int_y2;
    double m_Int_x2y;
    double m_Int_y2x;
    double m_Int_x3;
    double m_Int_y3;

  };


  class Planar2DPolygon {

   public:

    Planar2DPolygon();

    void Add(const Vector2d &point);

    double GetSurfaceIntegral(PolygonSurfaceIntegrals::IntegrandType type) const;

    void Update() const;

   private:
    void UpdateSurfaceIntegrals() const;

   private:
    std::vector<Vector2d> m_points;

    // Cache properties
    mutable bool c_is_up_to_date;
    mutable PolygonSurfaceIntegrals c_surface_integrals;

  };


// TODO: utiliser aabb d'eigen !

//  /**
//   * Class for a 2D Axis Aligned Bounding Box (AABB)
//   */
//  struct AABB2D { // TODO: placer cette classe ailleurs
//    double xmin = 0.;
//    double xmax = 0.;
//    double ymin = 0.;
//    double ymax = 0.;
//  };


  /**
     * Class for a polygon, based on a set of vertices.
     * Surface integrals over the surface delimited by the polygon can be computed
     */
  class Planar3DPolygon {

   public:

    explicit Planar3DPolygon(std::shared_ptr<Plane> plane);

    void AddPoint(const Vector3d &point);

    /// Get the plane passing through all the vertices, (be careful to check that your polygon is planar)
    /// \return plane related to the polygon
    std::shared_ptr<Plane> GetPlane() const;

//    /// Get the 2D Axis Aligned Bounding Box of the polygon
//    /// \return
//    AABB2D GetBoundingBox() const;

    Planar2DPolygon Get2DPolygon() const;

    bool IsEmpty() const {
      return m_points.empty();
    }

//   private:
    void Generate2DPolygon() const;

    void Update() const;

   private:
    std::shared_ptr<Plane> m_plane;
    std::vector<Vector3d> m_points;

    // Cache attributes
    mutable bool c_is_up_to_date;
    mutable Planar2DPolygon c_2D_polygon;

  };


  class Planar3DPolygonSet {

   public:
    explicit Planar3DPolygonSet(std::shared_ptr<Plane> common_plane) : m_common_plane(std::move(common_plane)) {}

    void AddPolygon(const Planar3DPolygon& polygon) {
      if (polygon.GetPlane() != m_common_plane) {
        std::cerr << "In a Planar3DPolygonSet every added polygon MUST share the same Plane" << std::endl;
        exit(EXIT_FAILURE);
      }
      m_poygons.push_back(polygon);
    }

    std::shared_ptr<Plane> GetPlane() const {
      return m_common_plane;
    }

   private:
    std::shared_ptr<Plane> m_common_plane;
    std::vector<Planar3DPolygon> m_poygons;

  };

}  // end namespace meshoui



#endif //MESHOUI_POLYGON_H

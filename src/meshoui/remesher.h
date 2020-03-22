//
// Created by frongere on 18/03/2020.
//

#ifndef MESHOUI_REMESHER_H
#define MESHOUI_REMESHER_H

#include <libmmgs.h>

#include "mesh.h"

// Temporaire
#include "meshoui/maths.h"
#include "vtkmesh.h"


#define MAX0(a, b)     (((a) > (b)) ? (a) : (b))
#define MAX3(a, b, c)  (((MAX0(a,b)) > c) ? (MAX0(a,b)) : c)


namespace meshoui {

//  namespace internal {
//
//    void OpenMeshToMMG5_pmesh(Mesh &omesh, MMG5_pMesh &mmg_mesh,
//                              MMG5_pSol mmg_sol) { // TODO: voir si on passe pas mmg_mesh et mmg_sol par ref...
//
//      // Setting the size of the mesh
//      MMGS_Set_meshSize(mmg_mesh, omesh.n_vertices(), omesh.n_faces(), 0); // C'est quoi na ???
//
//      int idx = 0;
//      auto vh_iter = omesh.vertices_begin();
//      for (; vh_iter != omesh.vertices_end(); ++vh_iter) {
//        idx++;
//        auto p = omesh.point(*vh_iter);
//        if (MMGS_Set_vertex(mmg_mesh, p[0], p[1], p[2], 0, idx) != 1) { // TODO: C'est quoi ref ?
//          std::cerr << "Unable to add vertex" << std::endl;
//          exit(EXIT_FAILURE);
//        }
//      }
//
//      std::array<int, 3> fidx = {0, 0, 0};
//      idx = 0;
//      auto fh_iter = omesh.faces_begin();
//      for (; fh_iter != omesh.faces_end(); ++fh_iter) {
//        idx++;
//
//        auto fv_iter = omesh.fv_ccwiter(*fh_iter);
//        int i = 0;
//        for (; fv_iter.is_valid(); ++fv_iter) {
////          if ((*fv_iter).idx() == 0) {
//////            std::cout << "0 iondexing" << std::endl;
//////          }
//          fidx[i] = (*fv_iter).idx() + 1; // FIXME: voir l'indexing de mmg (commence a 1 apparemment...)
//          i++;
//        }
//
//
//        if (MMGS_Set_triangle(mmg_mesh, fidx[0], fidx[1], fidx[2], 0, idx) != 1) { // TODO: C'est quoi ref ?
//          std::cerr << "Unable to add triangle" << std::endl;
//          exit(EXIT_FAILURE);
//        }
//      }
//
//    }
//  }


  class Remesher {

   public:
    Remesher();

    void SetHausdorffParam(const double& val);

    void SetMaxEdgeSize(const double& val);

    void SetConstantEdgeSize(const double& val);

    void SetAngleDetectionThreshold(const double& val_deg);

    int Apply(Mesh *mesh);

   private:

    using MMGMesh = MMG5_pMesh;
    using MMGSol = MMG5_pSol;

    void Initialize(MMGMesh& mmg_mesh, MMGSol& mmg_sol);

    void Finalize(MMGMesh& mmg_mesh, MMGSol& mmg_sol);

    static void OMesh2MMGMesh(Mesh &omesh, MMGMesh &mmg_mesh, MMGSol &mmg_sol);

    void ApplyParameters(MMGMesh& mmg_mesh, MMGSol& mmg_sol);

    static int Remesh(MMGMesh& mmg_mesh, MMGSol& mmg_sol);

    static void Check(MMGMesh& mmg_mesh, MMGSol& mmg_sol);

    static void ExtractSolution(Mesh &omesh, MMGMesh &mmg_mesh, MMGSol &mmg_sol);


   private:
    double m_hmax;
    double m_hausdorff_param;
    double m_constant_edge_size;
    double m_detection_angle;

  };

}  // end namespace meshoui



#endif //MESHOUI_REMESHER_H

//
// Created by frongere on 18/03/2020.
//

#include "remesher.h"

namespace meshoui {

  Remesher::Remesher() :
      m_constant_edge_size(5.),
      m_hausdorff_param(1.),
      m_hmax(20.),
      m_detection_angle(25) {}

  int Remesher::Apply(Mesh *mesh) {

    MMGMesh mmg_mesh = nullptr;
    MMGSol mmg_sol = nullptr;

    // mmgs Initialization
    Initialize(mmg_mesh, mmg_sol);

    // Generate a MMGMesh from meshoui::Mesh
    OMesh2MMGMesh(*mesh, mmg_mesh, mmg_sol);

    // Check if the number of given entities match with mesh size
    Check(mmg_mesh, mmg_sol);

    // Applying current remeshing parameters
    ApplyParameters(mmg_mesh, mmg_sol);

    // Remeshing
    int err_code = Remesh(mmg_mesh, mmg_sol);

    if (err_code == MMG5_STRONGFAILURE) {
      fprintf(stdout, "BAD ENDING OF MMGSLIB: UNABLE TO SAVE MESH\n");
      return err_code;

    } else if (err_code == MMG5_LOWFAILURE) {
      fprintf(stdout, "BAD ENDING OF MMGSLIB\n");
    }

    /**
     * Get the remeshed data
     */

    // The new remeshed meshoui mesh
    Mesh output_mesh;

    vtkSmartPointer<vtkPolyData> polydata;



    // Adding properties to get insights into feature edges and valence of each vertex
    auto feature_vertex_prop = output_mesh.CreateVertexProperty<bool>("is_feature_vertex");


    int nb_vertices, nb_triangles, nb_edges; // TODO: voir ce que nb_edges...

    if (MMGS_Get_meshSize(mmg_mesh, &nb_vertices, &nb_triangles, &nb_edges) != 1) exit(EXIT_FAILURE);


    /**
     * Vocabulaire:
     *
     * Un vertex is_required est un vertex qui doit etre conserve
     *
     * corners est un vertex qui
     *
     */

    /* Table to know if a vertex is corners */
    int *corners = (int *) calloc(nb_vertices + 1, sizeof(int));
    if (!corners) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    /* Table to know if a vertex/tetra/tria/edge is is_required */
    int *is_required = (int *) calloc(MAX3(nb_vertices, nb_triangles, nb_edges) + 1, sizeof(int));
    if (!is_required) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }
    /* Table to know if a coponant is corners and/or is_required */
    int *ridge_edges = (int *) calloc(nb_edges + 1, sizeof(int));
    if (!ridge_edges) {
      perror("  ## Memory problem: calloc");
      exit(EXIT_FAILURE);
    }


    /**
     * Get vertices
     */
    int nb_corners = 0;
    int nb_required = 0;
    for (int k = 1; k <= nb_vertices; k++) {
      Vector3d Point;
      int ref;

      // Note that each call to the following function increment an internal counter inside the MMGMesh (awful !)
      if (MMGS_Get_vertex(mmg_mesh, &(Point[0]), &(Point[1]), &(Point[2]),
                          &ref, &(corners[k]), &(is_required[k])) != 1) {
        exit(EXIT_FAILURE);
      }

//      // TODO: ici, on voudrait qualifier le vertex
//      if (mmg_mesh->point[k].tag & MG_CRN & corners[k]) {
//        std::cout << k << " is a corner" << std::endl;
//      }




      // add vertices to the output_mesh
      output_mesh.add_vertex(Point);
      if (corners[k]) nb_corners++;
      if (is_required[k]) nb_required++;
    }

    /**
     * Get Corners
     */
//      fprintf(inm,"\nCorners\n%d\n",nb_corners);
    for (int k = 1; k <= nb_vertices; k++) {
      if (corners[k]) {
//          fprintf(inm,"%d \n",k);
        // TODO: voir si on prend les corners
      }
    }

    /**
     * Required vertices
     */
//      fprintf(inm,"\nRequiredVertices\n%d\n",nb_required);
    for (int k = 1; k <= nb_vertices; k++) {
      if (is_required[k] != corners[k]) {
//        std::cerr << "corner which is not required: " << k << std::endl;
//          fprintf(inm,"%d \n",k);
        // TODO: voir si on traite cela
      }
    }
    free(corners);
    corners = NULL;

    /**
     * Get the triangles
     */
    nb_required = 0;
//      fprintf(inm,"\nTriangles\n%d\n",nb_triangles);
    for (int k = 1; k <= nb_triangles; k++) {
      int Tria[3];
      int ref;

      if (MMGS_Get_triangle(mmg_mesh, &(Tria[0]), &(Tria[1]), &(Tria[2]),
                            &ref, &(is_required[k])) != 1) { // TODO: Un triangle est required si ses 3 sommets le sont
        exit(EXIT_FAILURE);
      }

      // Add face to the output_mesh
      std::vector<Mesh::VertexHandle> face_vhandles;
      face_vhandles.clear();

      for (const auto &idx: Tria) {
        face_vhandles.emplace_back(Mesh::VertexHandle(idx - 1));
      }
      output_mesh.add_face(face_vhandles);

//        fprintf(inm,"%d %d %d %d \n",Tria[0],Tria[1],Tria[2],ref);
      if (is_required[k]) nb_required++;
    }

    /**
     * Required triangles
     */
//      fprintf(inm,"\nRequiredTriangles\n%d\n",nb_required);
    for (int k = 1; k <= nb_triangles; k++) {
      if (is_required[k]) {
//          fprintf(inm,"%d \n",k);
        // TODO: on fait quelque chose avec ca ?
      }
    }

    /**
     * Edges FIXME: c'est quoi ici ? feature edges ?
     */
    nb_required = 0;
    int nr = 0;
//      fprintf(inm,"\nEdges\n%d\n",nb_edges);
    for (int k = 1; k <= nb_edges; k++) {
      int Edge[2];
      int ref;

      if (MMGS_Get_edge(mmg_mesh, &(Edge[0]), &(Edge[1]), &ref,
                        &(ridge_edges[k]), &(is_required[k])) != 1) {
        exit(EXIT_FAILURE);
      }

//        fprintf(inm,"%d %d %d \n",Edge[0],Edge[1],ref);
      if (ridge_edges[k]) nr++;
      if (is_required[k]) nb_required++;
    }

    /**
     * Required edges
     */
//      fprintf(inm,"\nRequiredEdges\n%d\n",nb_required);
    for (int k = 1; k <= nb_edges; k++) {
      if (is_required[k]) {
//          fprintf(inm,"%d \n",k);
        // TODO: on fait quoi ?
      }
    }

    /**
     * Ridges
     */
//      fprintf(inm,"\nRidges\n%d\n",nr);
    for (int k = 1; k <= nb_edges; k++) {
      if (ridge_edges[k]) {
//          fprintf(inm,"%d \n",k);
        // TODO: On fait quoi ?
      }
    }

//      fprintf(inm,"\nEnd\n");
//      fclose(inm);

    free(is_required);
    is_required = nullptr;
    free(ridge_edges);
    ridge_edges = nullptr;

    // Freeing memory for mmgs
    Finalize(mmg_mesh, mmg_sol);

    // Replacing old mesh by remeshed one
    *mesh = output_mesh;


    return err_code;

  }

  void Remesher::Initialize(Remesher::MMGMesh &mmg_mesh, Remesher::MMGSol &mmg_sol) {
    MMGS_Init_mesh(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,
                   &mmg_mesh,
                   MMG5_ARG_ppMet,
                   &mmg_sol,
                   MMG5_ARG_end);
  }

  void Remesher::Finalize(Remesher::MMGMesh &mmg_mesh, Remesher::MMGSol &mmg_sol) {
    MMGS_Free_all(MMG5_ARG_start,
                  MMG5_ARG_ppMesh, &mmg_mesh, MMG5_ARG_ppMet, &mmg_sol,
                  MMG5_ARG_end);
  }

  void Remesher::SetHausdorffParam(const double &val) {
    m_hausdorff_param = val;
  }

  void Remesher::SetMaxEdgeSize(const double &val) {
    m_hmax = val;
  }

  void Remesher::SetConstantEdgeSize(const double &val) {
    m_constant_edge_size = val;
  }

  void Remesher::SetAngleDetectionThreshold(const double &val_deg) {
    m_detection_angle = val_deg;
  }

  void Remesher::OMesh2MMGMesh(Mesh &omesh, Remesher::MMGMesh &mmg_mesh, Remesher::MMGSol &mmg_sol) {
    // Setting the size of the mesh
    MMGS_Set_meshSize(mmg_mesh, omesh.n_vertices(), omesh.n_faces(), 0);

    int idx = 0;
    auto vh_iter = omesh.vertices_begin();
    for (; vh_iter != omesh.vertices_end(); ++vh_iter) {
      idx++;
      auto p = omesh.point(*vh_iter);
      if (MMGS_Set_vertex(mmg_mesh, p[0], p[1], p[2], 0, idx) != 1) { // TODO: C'est quoi ref ?
        std::cerr << "Unable to add vertex" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    std::array<int, 3> fidx = {0, 0, 0};
    idx = 0;
    auto fh_iter = omesh.faces_begin();
    for (; fh_iter != omesh.faces_end(); ++fh_iter) {
      idx++;

      auto fv_iter = omesh.fv_ccwiter(*fh_iter);
      int i = 0;
      for (; fv_iter.is_valid(); ++fv_iter) {
        fidx[i] = (*fv_iter).idx() + 1; // FIXME: voir l'indexing de mmg (commence a 1 apparemment...)
        i++;
      }

      if (MMGS_Set_triangle(mmg_mesh, fidx[0], fidx[1], fidx[2], 0, idx) != 1) { // TODO: C'est quoi ref ?
        std::cerr << "Unable to add triangle" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  void Remesher::ApplyParameters(MMGMesh &mmg_mesh, MMGSol &mmg_sol) {
    MMGS_Set_dparameter(mmg_mesh, mmg_sol, MMGS_DPARAM_hmax, m_hmax);
    MMGS_Set_dparameter(mmg_mesh, mmg_sol, MMGS_DPARAM_hausd, m_hausdorff_param);
    MMGS_Set_dparameter(mmg_mesh, mmg_sol, MMGS_DPARAM_hsiz, m_constant_edge_size);
    MMGS_Set_dparameter(mmg_mesh, mmg_sol, MMGS_DPARAM_angleDetection, m_detection_angle);
  }

  int Remesher::Remesh(Remesher::MMGMesh &mmg_mesh, Remesher::MMGSol &mmg_sol) {
    return MMGS_mmgslib(mmg_mesh, mmg_sol);
  }

  void Remesher::Check(Remesher::MMGMesh &mmg_mesh, Remesher::MMGSol &mmg_sol) {
    if (MMGS_Chk_meshData(mmg_mesh, mmg_sol) != 1) {
      std::cerr << "Check data failed" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  void Remesher::ExtractSolution(Mesh &omesh, Remesher::MMGMesh &mmg_mesh, Remesher::MMGSol &mmg_sol) {

  }

}  // end namespace meshoui

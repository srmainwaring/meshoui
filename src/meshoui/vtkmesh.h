// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_VTKMESH_H
#define MESHOUI_VTKMESH_H

#include "maths.h"
#include "mesh.h"

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCornerAnnotation.h>
#include <vtkTextProperty.h>
#include <vtkPolyDataMapper.h>
#include "vtkProperty.h"

namespace meshoui {

/**
  * Class for handling meshes using VTK.
  */
  class VTKMesh {

   public:

    enum WHERE {
      CELL,
      VERTEX
    };

    /// Constructor of the class.
    explicit VTKMesh(meshoui::Mesh &mesh);

    /// This function adds a field made of double to the polydata.
    void AddField(const std::string &name, const std::vector<double> &data, WHERE where);

    /// This function adds a field made of Vector 3d to the polydata.
    void AddField(const std::string &name, const std::vector<Vector3d> &data, WHERE where);

    /// This function writes the ouput *.vtp file.
    void Write(const std::string &meshfile);

    /// This function displays the mesh in a window.
    void Visualize();

    /// This function adds data of type T to either the faces or the vertices of the vtk mesh.
    template<typename T>
    void AddDynamicField(Mesh *mesh, const char *PropertyName, WHERE where) {

      std::vector<T> data;

      if (where == CELL) {
        auto property = OpenMesh::getProperty<OpenMesh::FaceHandle, T>(*mesh, PropertyName);
        data.reserve(mesh->n_faces());
        auto f_iter = mesh->faces_begin();
        for (; f_iter != mesh->faces_end(); ++f_iter) {
          data.push_back(mesh->property(*property, *f_iter));
        }
      } else if (where == VERTEX) {
        auto property = OpenMesh::getProperty<OpenMesh::VertexHandle, T>(*mesh, PropertyName);
        data.reserve(mesh->n_vertices());
        auto v_iter = mesh->vertices_begin();
        for (; v_iter != mesh->vertices_end(); ++v_iter) {
          data.push_back(mesh->property(*property, *v_iter));
        }
      }

      AddField(PropertyName, data, where);

    }

    /// This function adds the field of the face normals to the vtk mesh.
    void AddFaceNormalField(meshoui::Mesh *mesh);

    /// This function adds the field of the vertex normals to the vtk mesh.
    void AddVertexNormalField(meshoui::Mesh *mesh);

    /// This function adds the field of the edge lengths to the vtk mesh.
    void AddEdgeLengthField(meshoui::Mesh *mesh);

   private:

    /// This function fills the polydata.
    void FillPolyData(meshoui::Mesh &mesh);

   private:
    vtkSmartPointer<vtkPolyData> m_polydata; // Polydata structure of VTK.

  };

  void Show(Mesh &mesh);

  void Write_VTK(Mesh &mesh, const std::string &vtp_filename);


}  // end namespace meshoui

#endif // MESHOUI_VTKMESH_H

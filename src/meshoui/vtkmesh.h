// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#ifndef MESHOUI_VTKMESH_H
#define MESHOUI_VTKMESH_H

#include "math.h"
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

   private:

    /// This function fills the polydata.
    void FillPolyData(meshoui::Mesh &mesh);

   private:
    vtkSmartPointer<vtkPolyData> m_polydata; // Polydata structure of VTK.

  };

}  // end namespace meshoui

#endif // MESHOUI_VTKMESH_H

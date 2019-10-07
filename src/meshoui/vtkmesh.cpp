// ==========================================================================
// MeshOui
//
// Copyright (c) D-ICE Engineering.
// All rights reserved.
// ==========================================================================

#include "vtkmesh.h"


namespace meshoui {


  VTKMesh::VTKMesh(meshoui::Mesh &mesh) : m_polydata(vtkSmartPointer<vtkPolyData>::New()) {

    // Constructor of the class.

    // Creation of m_polydata.
    FillPolyData(mesh);

  }

  void VTKMesh::FillPolyData(meshoui::Mesh &mesh) {

    // This function fills the polydata.

    // Creation of the set of nodes.
    auto points = vtkSmartPointer<vtkPoints>::New();

    // Loop over the nodes.
    auto v_it = mesh.vertices_begin();
    for (; v_it != mesh.vertices_end(); ++v_it) {

      // Conversion from OpenMesh::point to Vector3d.
      Vector3d vertex = OpenMeshPointToVector3d<Vector3d>(mesh.point(*v_it));

      // Add a point to the structure the array of points.
      points->InsertNextPoint(vertex(0), vertex(1), vertex(2));

    }

    // Creation of the set of triangles.
    auto triangles = vtkSmartPointer<vtkCellArray>::New();

    // Create a triangle.
    auto triangle = vtkSmartPointer<vtkTriangle>::New();

    // Loop over the panels.
    int id;
    auto f_it = mesh.faces_begin();
    for (; f_it != mesh.faces_end(); ++f_it) {

      // Loop over the vertices of the triangle.
      id = 0;
      auto fv_it = mesh.fv_iter(*f_it);
      for (; fv_it.is_valid(); fv_it++) {
        triangle->GetPointIds()->SetId(id, (*fv_it).idx());
        id++;
      }

      // Check there are only three vertices per triangle.
      assert(id == 3);

      // Add the new panel.
      triangles->InsertNextCell(triangle);

    }

    // Add the nodes and the panels to polydata.
    m_polydata->SetPoints(points);
    m_polydata->SetPolys(triangles);

  }

  void VTKMesh::AddField(const std::string &name, const std::vector<double> &data, WHERE where) {

    // This function adds a field to the polydata.

    if (where == VERTEX) {
      assert(data.size() == m_polydata->GetNumberOfPoints());
    } else if (where == CELL) {
      assert(data.size() == m_polydata->GetNumberOfCells());
    }


    auto field = vtkSmartPointer<vtkDoubleArray>::New();
    field->SetNumberOfComponents(1);
    field->SetName(name.c_str());


    vtkSmartPointer<vtkUnsignedCharArray> colors =
        vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");


    for (const auto &val : data) {
      field->InsertNextValue(val);
      colors->InsertNextValue(val);
    }

    if (where == VERTEX) {
      m_polydata->GetPointData()->AddArray(field);

      m_polydata->GetPointData()->SetScalars(colors);

    } else if (where == CELL) {
      m_polydata->GetCellData()->AddArray(field);
    }

  }

  void VTKMesh::Write(const std::string &meshfile) {

    // This function writes the ouput *.vtp file.

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(meshfile.c_str());
    writer->SetInputData(m_polydata);
    writer->Write();
  }

  const vtkSmartPointer<vtkPolyData> &VTKMesh::polydata(){

    // Getter for m_polydata.

    return m_polydata;

  }

} // end namespace meshoui

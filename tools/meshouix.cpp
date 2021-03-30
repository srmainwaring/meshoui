//
// Created by frongere on 25/09/2020.
//

#include <meshoui/meshoui.h>
#include <cxxopts.hpp>

#include <filesystem>

namespace fs = std::filesystem;

bool has_obj_extension(const std::string &file) {
  return fs::path(file).extension() == ".obj";
}

int main(int argc, char *argv[]) {

  cxxopts::Options options("meshouix",
                           "A D-ICE program to apply remeshing using mmgs library");

  options.add_options()
      ("i,input", "Input mesh file with obj file format", cxxopts::value<std::string>())
      ("o,output", "Ouput mesh file with remeshing applied with obj file format", cxxopts::value<std::string>())

      ("c,clip", "First clip the input mesh by the plane Oxy", cxxopts::value<bool>()->default_value("false"))


      ("a,angle", "Angle detection threshold in degrees", cxxopts::value<double>()->default_value("40"))
      ("p,haussdorf", "Haussdorf Parameter", cxxopts::value<double>()->default_value("0.5"))
      ("e,edge", "Target edge size", cxxopts::value<double>()->default_value("1"))

      ("n,no-writing", "No output file writing", cxxopts::value<bool>()->default_value("false"))

      ("show", "Show result in vtk figure", cxxopts::value<bool>()->default_value("false"))

      ("h,help", "Print help");

  // TODO: ajouter une option de clipping par le plan horizontal en z=0

  options.parse_positional({"input", "output"});

  // Parsing the command line arguments
  auto result = options.parse(argc, argv);

  std::cout << argv[0] << std::endl;
  std::cout << "=======================================================" << std::endl;
  std::cout << "MeshouiX quick remeshing tool" << std::endl;
  std::cout << "Property of D-ICE engineering" << std::endl;
  std::cout << "=======================================================" << std::endl;
  std::cout << std::endl;

  // If help options, print it and exit
  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    exit(EXIT_SUCCESS);
  }

  // Getting options values
  auto angle_threshold = result["angle"].as<double>();
  auto haussdorf_param = result["haussdorf"].as<double>();
  auto edge_length = result["edge"].as<double>();
  auto show = result["show"].as<bool>();

  // Input file as a mandatory positional argument
  std::string input_file;
  if (result.count("input")) {
    input_file = result["input"].as<std::string>();
    if (!has_obj_extension(input_file)) {
      std::cerr << "Input file (" << input_file << ") must have a .obj extension" << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cerr << "First command line argument must be a mesh file with obj file format" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "Reading file " << input_file << std::endl;
  meshoui::Mesh mesh(input_file);

  if (result.count("clip")) {
    std::cout << "Clipping the mesh by plane Oxy" << std::endl;
    auto plane = std::make_shared<meshoui::Plane>(meshoui::Vector3d{0., 0., 0.},
                                                  meshoui::Vector3d{0., 0.,
                                                                    1.}); // FIXME: pourquoi ce constructeur ne fonctionen pas ???
    meshoui::Clipper<meshoui::ClippingPlane> clipper(std::make_shared<meshoui::ClippingPlane>(plane));
    clipper.ClipIt(mesh);
  }


  // Output for current options
  std::cout << std::endl;
  std::cout << "Remeshing with the following parameters:" << std::endl;
  std::cout << "\t* Angle detection threshold (deg)    " << angle_threshold << std::endl;
  std::cout << "\t* Haussdorf parameter                " << haussdorf_param << std::endl;
  std::cout << "\t* Constant edge length (m)           " << edge_length << std::endl;
  std::cout << std::endl;

  // Instance of the remesher
  meshoui::Remesher remesher;
  remesher.SetAngleDetectionThreshold(angle_threshold);
  remesher.SetHausdorffParam(haussdorf_param);
  remesher.SetConstantEdgeSize(edge_length);

  remesher.RemeshIt(mesh);


  // Displaying the remeshed result into a VTK window
  if (show) {
    meshoui::Show(mesh);
  }


  if (!result.count("no-writing")) {

    std::string output_file;
    if (result.count("output")) {
      output_file = result["output"].as<std::string>();
      if (!has_obj_extension(input_file)) {
        std::cerr << "Output file (" << output_file << ") must have a .obj extension" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      // No output file name has been given, building one from the input file
      output_file = (std::string) fs::path(input_file).stem() + "_rem.obj";
    }

    // TODO: donner la possibilite de ne rien ecrire... (-nw)
    std::cout << "Writing remeshed file to " << output_file << std::endl;
    meshoui::Write_OBJ(mesh, output_file);
  }


  return 0;
}
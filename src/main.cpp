#include "geodesicdistance.hpp"
#include "vcgtrimesh.hpp"
#include <filesystem>
#include <igl/avg_edge_length.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/heat_geodesics.h>
#include <igl/isolines_map.h>
#include <igl/opengl/create_shader_program.h>
#include <igl/opengl/destroy_shader_program.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangulated_grid.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>
void set_colormap(igl::opengl::glfw::Viewer &viewer) {
  const int num_intervals = 30;
  Eigen::MatrixXd CM(num_intervals, 3);
  // Colormap texture
  for (int i = 0; i < num_intervals; i++) {
    double t = double(num_intervals - i - 1) / double(num_intervals - 1);
    CM(i, 0) = std::max(std::min(2.0 * t - 0.0, 1.0), 0.0);
    CM(i, 1) = std::max(std::min(2.0 * t - 1.0, 1.0), 0.0);
    CM(i, 2) = std::max(std::min(6.0 * t - 5.0, 1.0), 0.0);
  }
  igl::isolines_map(Eigen::MatrixXd(CM), CM);
  viewer.data().set_colormap(CM);
}

int main(int argc, char *argv[]) {
  // Create the peak height field
  Eigen::MatrixXi F;
  Eigen::MatrixXd V;
  // Performance table
  //  const std::string meshFilepath = "/home/iason/Models/vaselion.ply";
  //  const std::string meshFilepath = "/home/iason/Models/bimba.ply";
  //  const std::string meshFilepath = "/home/iason/Models/buste.ply";
  //  const std::string meshFilepath = "/home/iason/Models/fertility.ply";
  //  const std::string meshFilepath = "/home/iason/Models/dragon.ply";
  //  const std::string meshFilepath =
  //  "/home/iason/Models/fertility_highRes.obj"; const std::string meshFilepath
  //  = "/home/iason/Models/horse_1M.ply";

  const std::string meshFilepath = "/home/iason/Models/Armadillo_drilled.obj";

  //  const std::string meshFilepath = "/home/iason/Models/bunny.obj";
  //  const std::string meshFilepath = "/home/iason/Models/bunny_refined.obj";

  //  const std::string meshFilepath = "/home/iason/Models/torus.ply";

  //  const std::string meshFilepath = "/home/iason/Models/bunny.ply";
  //  const std::string meshFilepath =
  //  "/home/iason/Models/bunny_remeshed_300k.ply"	;
  //  const std::string meshFilepath = "/home/iason/Models/Armadillo.ply";

  //  const std::string meshFilepath = "/home/iason/Models/bunny.obj";
  //  const std::string meshFilepath =
  //  "/home/iason/Models/bunny_remeshed_300k.ply";
  //  const std::string meshFilepath = "/home/iason/Models/fertility_20k.ply";
  //    const std::string meshFilepath = "/home/iason/Models/bimba.ply";

  if (!std::filesystem::exists(std::filesystem::path(meshFilepath))) {
    std::cerr << meshFilepath << " does not exist." << std::endl;
    std::terminate();
  }
  igl::read_triangle_mesh(meshFilepath, V, F);

  // Precomputation
  igl::HeatGeodesicsData<double> data;
  double hSquared = std::pow(igl::avg_edge_length(V, F), 2);
  double m = 1;
  double t = m * hSquared;
  const auto precompute = [&]() {
    if (!igl::heat_geodesics_precompute(V, F, t, data)) {
      std::cerr << "Error: heat_geodesics_precompute failed." << std::endl;
      exit(EXIT_FAILURE);
    };
  };
  //  precompute();

  const double bboxDiag = igl::bounding_box_diagonal(V);
  VCGTriMesh mesh(meshFilepath);
  vcg::tri::UpdatePosition<VCGTriMesh>::Scale(mesh, 100 / bboxDiag);
  GeodesicDistance geodesicDistanceComputer(mesh);
  std::unordered_set<VertexIndex> sourcesVi;

  igl::opengl::glfw::Viewer viewer;
  bool useIgl = false;
  bool down_on_mesh = false;
  const auto update = [&]() -> bool {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                 viewer.core().proj, viewer.core().viewport, V,
                                 F, fid, bc)) {
      Eigen::VectorXd D;
      // if big mesh, just use closest vertex. Otherwise, blend distances to
      // vertices of face using barycentric coordinates.
      if (F.rows() > 100000) {
        // 3d position of hit
        const Eigen::RowVector3d m3 = V.row(F(fid, 0)) * bc(0) +
                                      V.row(F(fid, 1)) * bc(1) +
                                      V.row(F(fid, 2)) * bc(2);
        int cid = 0;
        Eigen::Vector3d((V.row(F(fid, 0)) - m3).squaredNorm(),
                        (V.row(F(fid, 1)) - m3).squaredNorm(),
                        (V.row(F(fid, 2)) - m3).squaredNorm())
            .minCoeff(&cid);
        const int vid = F(fid, cid);
        if (useIgl) {
          igl::heat_geodesics_solve(
              data, (Eigen::VectorXi(1, 1) << vid).finished(), D);
        } else {
          sourcesVi.clear();
          sourcesVi.insert(vid);
          std::unordered_map<VertexIndex, double> distanceMap;
          D = geodesicDistanceComputer.computeGeodesicDistances(sourcesVi,
                                                                distanceMap);
        }
      } else {
        D = Eigen::VectorXd::Zero(V.rows());
        for (int cid = 0; cid < 3; cid++) {
          const int vid = F(fid, cid);
          Eigen::VectorXd Dc;
          if (useIgl) {
            igl::heat_geodesics_solve(
                data, (Eigen::VectorXi(1, 1) << vid).finished(), Dc);
          } else {
            sourcesVi.clear();
            sourcesVi.insert(vid);
            std::unordered_map<VertexIndex, double> distanceMap;
            Dc = geodesicDistanceComputer.computeGeodesicDistances(sourcesVi,
                                                                   distanceMap);
          }
          D += Dc * bc(cid);
        }
      }
      viewer.data().set_data(D);
      return true;
    }
    return false;
  };
  viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer &viewer, int,
                                   int) -> bool {
    if (update()) {
      down_on_mesh = true;
      return true;
    }
    return false;
  };
  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &viewer, int,
                                   int) -> bool {
    if (down_on_mesh) {
      update();
      return true;
    }
    return false;
  };
  viewer.callback_mouse_up = [&down_on_mesh](igl::opengl::glfw::Viewer &viewer,
                                             int, int) -> bool {
    down_on_mesh = false;
    return false;
  };
  std::cout << R"(Usage:
    [click]  Click on shape to pick new geodesic distance source
    ,/.      Decrease/increase t by factor of 10.0
    D,d      Toggle using intrinsic Delaunay discrete differential operators
  )";

  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer & /*viewer*/,
                                    unsigned int key, int mod) -> bool {
    switch (key) {
    default:
      return false;
    case 'c':
    case 'C':
      useIgl = !useIgl;
      if (useIgl) {
        std::cout << "Using igl" << std::endl;
        precompute();
      } else {
        std::cout << "Using my implemenetation" << std::endl;
        geodesicDistanceComputer.setMFactor(m);
      }
      update();
      break;
    case 'D':
    case 'd':
      data.use_intrinsic_delaunay = !data.use_intrinsic_delaunay;
      std::cout << (data.use_intrinsic_delaunay ? "" : "not ")
                << "using intrinsic delaunay..." << std::endl;
      precompute();
      update();
      break;
    case '.':
    case ',':
      m *= (key == '.' ? 10.0 : 0.1);
      std::cout << "m: " << m << std::endl;
      t = m * hSquared;
      if (useIgl) {
        precompute();
      } else {
        geodesicDistanceComputer.setMFactor(m);
      }
      update();
      break;
    }
    return true;
  };

  // Show mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_data(Eigen::VectorXd::Zero(V.rows()));

  set_colormap(viewer);
  viewer.data().show_lines = false;
  viewer.launch();
}

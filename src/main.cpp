#include "geodesicdistance.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "vcgtrimesh.hpp"
#include <filesystem>

int main() {
  const std::string meshFilepath = "/home/iason/Models/buste.ply";
  //  const std::string meshFilepath = "/home/iason/Models/Greek_Sculpture.off";
  //  const std::string meshFilepath = "/home/iason/Models/fertility.ply";
  //  const std::string meshFilepath = "/home/iason/Models/bunny.obj";
  VCGTriMesh m(meshFilepath);

  polyscope::init();
  const std::string meshName = std::filesystem::path(meshFilepath).stem();
  polyscope::registerSurfaceMesh(meshName, m.getVertices(), m.getFaces());
  // Compute geodesic distances
  GeodesicDistance geodesicDistanceComputer(m);
  std::cout << "Built" << std::endl;
  std::unordered_map<VertexIndex, double> distanceMap;
  Eigen::VectorXd distances =
      geodesicDistanceComputer.computeGeodesicDistances({0}, distanceMap);
  std::cout << "Computed distances" << std::endl;

  polyscope::getSurfaceMesh(meshName)->addVertexDistanceQuantity(meshName,
                                                                 distances);

  polyscope::show();
}

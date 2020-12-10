#include "geodesicdistance.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "vcgtrimesh.hpp"
#include <filesystem>

int main() {
  //  const std::string meshFilepath = "/home/iason/Models/buste.ply";
  //  const std::string meshFilepath = "/home/iason/Models/Greek_Sculpture.off";
  //  const std::string meshFilepath = "/home/iason/Models/fertility.ply";
  //  const std::string meshFilepath = "/home/iason/Models/bunny_low.obj";
  const std::string meshFilepath = "/home/iason/Models/fertility_highRes.obj";
  VCGTriMesh m(meshFilepath);

  polyscope::init();
  const std::string meshName = std::filesystem::path(meshFilepath).stem();
  polyscope::registerSurfaceMesh(meshName, m.getVertices(), m.getFaces());
  // Compute geodesic distances
  GeodesicDistance geodesicDistanceComputer(m);
  const size_t viSource = 0;
  std::cout << "source: " << m.vert[viSource].cP()[0] << " "
            << m.vert[viSource].cP()[1] << " " << m.vert[viSource].cP()[2]
            << std::endl;
  std::unordered_map<VertexIndex, double> distanceMap;
  Eigen::VectorXd distances = geodesicDistanceComputer.computeGeodesicDistances(
      {viSource}, distanceMap);
  for (int i = 0; i < m.VN(); i++) {
    std::cout << "v" << i << "  is at distance " << distances(i)
              << " to v" + std::to_string(viSource) << std::endl;
  }
  std::cout << "Computed distances" << std::endl;

  polyscope::getSurfaceMesh(meshName)->addVertexDistanceQuantity(meshName,
                                                                 distances);

  polyscope::show();
}

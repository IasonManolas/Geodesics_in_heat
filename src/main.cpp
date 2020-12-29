#include "geodesicdistance.hpp"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "vcgtrimesh.hpp"
#include <filesystem>

int main() {
  //  const std::string meshFilepath = "/home/iason/Models/buste.ply";
  //  const std::string meshFilepath = "/home/iason/Models/Greek_Sculpture.off";
  const std::string meshFilepath = "/home/iason/Models/fertility.ply";
  //  const std::string meshFilepath = "/home/iason/Models/bunny_low.obj";
  //  const std::string meshFilepath =
  //  "/home/iason/Models/fertility_highRes.obj";
  VCGTriMesh m(meshFilepath);

  polyscope::init();
  const std::string meshName = std::filesystem::path(meshFilepath).stem();
  polyscope::registerSurfaceMesh(meshName, m.getVertices(), m.getFaces());
  // Compute geodesic distances
  GeodesicDistance geodesicDistanceComputer(m);
  const size_t viSource = 119612;
  std::cout << "source: " << m.vert[viSource].cP()[0] << " "
            << m.vert[viSource].cP()[1] << " " << m.vert[viSource].cP()[2]
            << std::endl;
  std::unordered_map<VertexIndex, double> distanceMap;
  std::unordered_set<VertexIndex> sourcesVi;
  sourcesVi.insert(viSource);
  Eigen::VectorXd distances =
      geodesicDistanceComputer.computeGeodesicDistances(sourcesVi, distanceMap);
  //  for (int i = 0; i < m.VN(); i++) {
  //    std::cout << "v" << i << "  is at distance " << distances(i)
  //              << " to v" + std::to_string(viSource) << std::endl;
  //  }
  std::cout << "Computed distances" << std::endl;

  polyscope::getSurfaceMesh(meshName)->addVertexDistanceQuantity(meshName,
                                                                 distances);
  std::vector<std::array<double, 3>> sourcePointCoords(sourcesVi.size());
  int sourceIndex = 0;
  for (auto sourceVi : sourcesVi) {
    sourcePointCoords[sourceIndex][0] = m.vert[sourceVi].cP()[0];
    sourcePointCoords[sourceIndex][1] = m.vert[sourceVi].cP()[1];
    sourcePointCoords[sourceIndex][2] = m.vert[sourceVi].cP()[2];
    sourceIndex++;
  }
  polyscope::registerPointCloud("Source points", sourcePointCoords);

  polyscope::show();
}

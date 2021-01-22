#include "geodesicdistance.hpp"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "vcgtrimesh.hpp"
#include <filesystem>

void addGaussianNoise(VCGTriMesh &m) {
  std::normal_distribution<double> dist(0, 1);
  std::default_random_engine generator;
  for (int vi = 0; vi < m.VN(); vi++) {
    m.vert[vi].P()[0] = m.vert[vi].P()[0] + dist(generator);
    m.vert[vi].P()[1] = m.vert[vi].P()[1] + dist(generator);
    m.vert[vi].P()[2] = m.vert[vi].P()[2] + dist(generator);
  }
}
int main() {
  //  const std::string meshFilepath = "/home/iason/Models/Armadillo.ply";
  //  const std::string meshFilepath =
  //  "/home/iason/Models/Armadillo_drilled.ply"; const std::string meshFilepath
  //  = "/home/iason/Models/buste.ply"; const std::string meshFilepath =
  //  "/home/iason/Models/buste_150K.ply"; const std::string meshFilepath =
  //  "/home/iason/Models/Greek_Sculpture.off"; const std::string meshFilepath =
  //  "/home/iason/Models/fertility.ply"; const std::string meshFilepath =
  //  "/home/iason/Models/bunny.ply"; const std::string meshFilepath =
  //  "/home/iason/Models/horse.ply"; const std::string meshFilepath =
  //  "/home/iason/Models/horse_1M.ply"; const std::string meshFilepath =
  //  "/home/iason/Models/horse_6M.ply"; const std::string meshFilepath =
  //  "/home/iason/Models/vaselion.ply";
  //    const std::string meshFilepath = "/home/iason/Models/bimba.ply";
  //  const std::string meshFilepath = "/home/iason/Models/bimba_100K.ply";
  //  const std::string meshFilepath = "/home/iason/Models/happy_remeshed.ply";
  //  const std::string meshFilepath = "/home/iason/Models/dragon.ply";
  //    const std::string meshFilepath = "/home/iason/Models/bunny_low.ply";
  const std::string meshFilepath = "/home/iason/Models/bunny_low_drilled.ply";

  //  const std::string meshFilepath = "/home/iason/Models/tetrahedron.ply";
  //  const std::string meshFilepath =
  //  "/home/iason/Models/fertility_highRes.obj";
  VCGTriMesh m(meshFilepath);
  const std::string meshName = std::filesystem::path(meshFilepath).stem();
  std::cout << "Number of verts:" << m.VN() << std::endl;
  std::cout << "Number of triangles:" << m.FN() << std::endl;

  // Compute geodesic distances
  auto begin = std::chrono::high_resolution_clock::now();
  GeodesicDistance geodesicDistanceComputer(m);
  auto end = std::chrono::high_resolution_clock::now();
  auto precomputeTime =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  std::cout << "Precompute in " << precomputeTime.count() / 1000.0 << " seconds"
            << std::endl;
  const size_t viSource = 0; // fert
  //  const size_t viSource = 119612; // fert
  //  const size_t viSource = 44023; // buste
  //  const size_t viSource = 13431; // horse
  //  const size_t viSource = 7174; // bimba
  //  const size_t viSource = 33295; // vaselion
  //  const size_t viSource = 287; // bunny
  //    const size_t viSource = 2176; // greek sculp
  //  const size_t viSource = 88459; // Armadillo
  //  std::cout << "source: " << m.vert[viSource].cP()[0] << " "
  //            << m.vert[viSource].cP()[1] << " " << m.vert[viSource].cP()[2]
  //            << std::endl;
  std::unordered_map<VertexIndex, double> distanceMap;
  std::unordered_set<VertexIndex> sourcesVi;
  sourcesVi.insert(viSource);
  Eigen::VectorXd distances =
      geodesicDistanceComputer.computeGeodesicDistances(sourcesVi, distanceMap);
  //  for (int i = 0; i < m.VN(); i++) {
  //    std::cout << "v" << i << "  is at distance " << distances(i)
  //              << " to v" + std::to_string(viSource) << std::endl;
  //  }
  end = std::chrono::high_resolution_clock::now();
  auto elapsedTotal =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
  std::cout << "Solving in "
            << elapsedTotal.count() / 1000.0 - precomputeTime.count() / 1000.0
            << " seconds" << std::endl;
  std::cout << "Computed distances in " << elapsedTotal.count() / 1000.0
            << " seconds" << std::endl;
  polyscope::init();
  polyscope::registerSurfaceMesh(meshName, m.getVertices(), m.getFaces());

  polyscope::getSurfaceMesh(meshName)->addVertexDistanceQuantity(meshName,
                                                                 distances);

  polyscope::getSurfaceMesh(meshName)->addVertexScalarQuantity("Colors",
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

  VCGTriMesh m_noisy;
  vcg::tri::Append<VCGTriMesh, VCGTriMesh>::MeshCopy(m_noisy, m);
  addGaussianNoise(m_noisy);
  polyscope::registerSurfaceMesh(meshName + "_withNoise", m_noisy.getVertices(),
                                 m_noisy.getFaces());
  GeodesicDistance geodesicDistanceComputer_noise(m_noisy);
  Eigen::VectorXd distances_noise =
      geodesicDistanceComputer_noise.computeGeodesicDistances(sourcesVi,
                                                              distanceMap);
  polyscope::getSurfaceMesh(meshName)->addVertexDistanceQuantity(
      meshName, distances_noise);

  polyscope::show();
}

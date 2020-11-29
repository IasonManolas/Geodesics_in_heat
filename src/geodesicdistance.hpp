#ifndef GEODESICDISTANCE_HPP
#define GEODESICDISTANCE_HPP
#include "vcgtrimesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <unordered_set>

using ConstVCGTriMesh = VCGTriMesh;
class GeodesicDistance {
  ConstVCGTriMesh &m;
  double m_timestep;
  Eigen::SparseVector<double> m_kroneckerDelta;
  Eigen::VectorXd m_solved_u;
  std::vector<VCGTriMesh::CoordType> m_unitGradient;
  Eigen::VectorXd m_divergence;

  bool isInitialized{false};

  void init();

public:
  GeodesicDistance(ConstVCGTriMesh &m);
  void setMesh(const VCGTriMesh &m);

  Eigen::VectorXd
  computeGeodesicDistances(const std::unordered_set<VertexIndex> &sourcesVi,
                           std::unordered_map<VertexIndex, double> &distances);

private:
  void updateKroneckerDelta(const std::unordered_set<VertexIndex> &sourcesVi);
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> la;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> la_cotan;
  void solveCotanLaplace();
  void computeUnitGradient();
  void computeDivergence();
  Eigen::VectorXd solvePhi(const std::unordered_set<VertexIndex> &sourcesVi);
};

#endif // GEODESICDISTANCE_HPP

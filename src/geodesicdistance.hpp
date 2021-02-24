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

  void precompute();

public:
  GeodesicDistance(ConstVCGTriMesh &m);
  void setMesh(ConstVCGTriMesh &m);

  Eigen::VectorXd
  computeGeodesicDistances(const std::unordered_set<VertexIndex> &sourcesVi,
                           std::unordered_map<VertexIndex, double> &distances);

  void setMFactor(double m);

  void setMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

private:
  void updateKroneckerDelta(const std::unordered_set<VertexIndex> &sourcesVi);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> la;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> la_cotan;
  void solveCotanLaplace();
  void computeUnitGradient();
  void computeDivergence();
  Eigen::VectorXd solvePhi(const std::unordered_set<VertexIndex> &sourcesVi);
  Eigen::SparseMatrix<double> Lc;
  Eigen::SparseMatrix<double> A;

  double averageEdgeLength;
  int debugVi{8};
  void init();
};

#endif // GEODESICDISTANCE_HPP

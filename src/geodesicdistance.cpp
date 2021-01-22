#include "geodesicdistance.hpp"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "vcg/complex/algorithms/mesh_to_matrix.h"
#include "vcg/complex/algorithms/stat.h"
#include <unordered_set>

GeodesicDistance::GeodesicDistance(ConstVCGTriMesh &m) : m(m) {
  vcg::tri::RequirePerFaceNormal(m);
  vcg::tri::UpdateNormal<VCGTriMesh>::PerFaceNormalized(m);
  init();
}

void GeodesicDistance::init() {
  assert(m.VN() != 0 && m.FN() != 0);
  std::cout << "Initializing.." << std::endl;

  // compute cotangent operator L_C
  std::cout << "Computing cotangent matrix.." << std::endl;
  std::vector<std::pair<int, int>> laplacianMatrixIndices;
  std::vector<double> laplacianMatrixEntries;
  vcg::tri::MeshToMatrix<VCGTriMesh>::GetLaplacianMatrix(
      m, laplacianMatrixIndices, laplacianMatrixEntries, true, 1, false);
  assert(laplacianMatrixIndices.size() == laplacianMatrixEntries.size());
  std::cout << "Populating sparse cotangent matrix.." << std::endl;
  Eigen::SparseMatrix<double> Lc(m.VN(), m.VN());
  const size_t numberOfNonZeroEntriesLc = laplacianMatrixEntries.size();
  std::vector<Eigen::Triplet<double>> LcTriplets(numberOfNonZeroEntriesLc);
  for (size_t nonZeroEntryIndex = 0;
       nonZeroEntryIndex < numberOfNonZeroEntriesLc; nonZeroEntryIndex++) {
    LcTriplets[nonZeroEntryIndex] =
        Eigen::Triplet<double>(laplacianMatrixIndices[nonZeroEntryIndex].first,
                               laplacianMatrixIndices[nonZeroEntryIndex].second,
                               laplacianMatrixEntries[nonZeroEntryIndex] * 0.5);
  }
  Lc.setFromTriplets(LcTriplets.begin(), LcTriplets.end());
  //  Eigen::MatrixXd LcDense(Lc);

  // compute mass matrix A
  std::cout << "Computing mass matrix.." << std::endl;
  std::vector<std::pair<int, int>> massMatrixIndices;
  std::vector<double> massMatrixEntries;
  vcg::tri::MeshToMatrix<VCGTriMesh>::MassMatrixEntry(m, massMatrixIndices,
                                                      massMatrixEntries, false);
  assert(massMatrixIndices.size() == massMatrixEntries.size());
  Eigen::VectorXd A_vec(m.VN());
  Eigen::SparseMatrix<double> A(m.VN(), m.VN());
  for (size_t vi = 0; vi < m.VN(); vi++) {
    assert(massMatrixIndices[vi].first == vi &&
           massMatrixIndices[vi].second == vi);
    A_vec(vi) = massMatrixEntries[vi];
  }

  A = A_vec.asDiagonal();

  // compute timestep
  const double averageEdgeLength =
      vcg::tri::Stat<VCGTriMesh>::ComputeEdgeLengthSum(m) / m.EN();
  m_timestep = averageEdgeLength * averageEdgeLength;

  // factor first equation
  std::cout << "Prefactoring the heat equation.." << std::endl;
  Eigen::SparseMatrix<double> B, B0;
  B0 = m_timestep * Lc;
  B = A /*..toDenseMatrix()*/ + B0; // shouldnt this be A-B0?
                                    //  Eigen::MatrixXd BDense(B);
  la.compute(B);
  if (la.info() != Eigen::Success) {
    std::cerr << "Prefactoring the heat equation failed." << std::endl;
    std::terminate();
  }

  // factor the poisson problem
  std::cout << "Prefactoring the poisson's equation.." << std::endl;
  la_cotan.compute(Lc);
  if (la_cotan.info() != Eigen::Success) {
    std::cerr << "Prefactoring the poisson equation failed." << std::endl;
    std::terminate();
  }
  isInitialized = true;
}

void GeodesicDistance::updateKroneckerDelta(
    const std::unordered_set<VertexIndex> &sourcesVi) {
  Eigen::SparseVector<double> delta(m.VN());
  delta.reserve(sourcesVi.size());
  for (VertexIndex vi : sourcesVi) {
    delta.insert(vi) = 1;
  }

  m_kroneckerDelta = delta;
}

void GeodesicDistance::solveCotanLaplace() {
  m_solved_u = la.solve(m_kroneckerDelta);
}

void GeodesicDistance::computeUnitGradient() {
  if (m_unitGradient.empty()) {
    m_unitGradient.resize(m.FN());
  }
  for (size_t fi = 0; fi < m.FN(); fi++) {
    VCGTriMesh::FaceType &f = m.face[fi];
    const double faceArea = vcg::DoubleArea(m.face[fi]) / 2;
    const VCGTriMesh::VertexType &v0 = *f.cV(0);
    const VCGTriMesh::VertexType &v1 = *f.cV(1);
    const VCGTriMesh::VertexType &v2 = *f.cV(2);
    double u0 = m_solved_u(m.getIndex(v0));
    double u1 = m_solved_u(m.getIndex(v1));
    double u2 = m_solved_u(m.getIndex(v2));

    double r_mag = 1. / std::max(std::max(u0, u1), u2);
    if (std::isfinite(r_mag)) {
      u0 *= r_mag;
      u1 *= r_mag;
      u2 *= r_mag;
    }
    auto vij = v1.cP() - v0.cP();
    auto pj = v1.cP();
    auto pk = v2.cP();
    VCGTriMesh::CoordType faceNormal = f.cN();
    VCGTriMesh::CoordType sum = (faceNormal ^ (v1.cP() - v0.cP())) * u2;
    sum += (faceNormal ^ (v2.cP() - v1.cP())) * u0;
    sum += (faceNormal ^ (v0.cP() - v2.cP())) * u1;
    sum /= faceArea;
    m_unitGradient[fi] = sum.Normalize();
  }
}

void GeodesicDistance::computeDivergence() {
  if (m_divergence.rows() == 0) {
    m_divergence.resize(m.VN());
  }
  m_divergence << Eigen::VectorXd::Zero(m.VN());

  // open a file in write mode.
  std::ofstream outfile;
  //  outfile.open("divergence.dat");
  for (size_t fi = 0; fi < m.FN(); fi++) {
    const VCGTriMesh::FaceType &f = m.face[fi];
    const VCGTriMesh::VertexPointer &v0 = f.cV(0);
    const VCGTriMesh::VertexPointer &v1 = f.cV(1);
    const VCGTriMesh::VertexPointer &v2 = f.cV(2);
    VCGTriMesh::CoordType e01(v1->cP() - v0->cP());
    VCGTriMesh::CoordType e12(v2->cP() - v1->cP());
    VCGTriMesh::CoordType e20(v0->cP() - v2->cP());

    const double norm_cross = vcg::DoubleArea(m.face[fi]);

    const double cotan0 = (e01 * (-e20)) / norm_cross;
    const double cotan1 = (e12 * (-e01)) / norm_cross;
    const double cotan2 = (e20 * (-e12)) / norm_cross;

    const VCGTriMesh::CoordType &X_fi = m_unitGradient[fi];
    // continue here
    //    auto i_entry = ((X_fi * e01) * cotan2 + (X_fi * (-e20)) * cotan1);
    //    auto j_entry = ((X_fi * e12) * cotan0 + (X_fi * (-e01)) * cotan2);
    //    auto k_entry = ((X_fi * e20) * cotan1 + (X_fi * (-e12)) * cotan0);
    //    outfile << "face " << fi << std::endl;
    //    outfile << "i entry " << i_entry << std::endl;
    //    outfile << "j entry " << j_entry << std::endl;
    //    outfile << "k entry " << k_entry << std::endl;
    //    auto vi0 = m.getIndex(v0);
    //    auto vi1 = m.getIndex(v1);
    //    auto vi2 = m.getIndex(v2);
    //    if (vi0 == 8) {
    //      std::cout << "fi " << fi << " " << i_entry << std::endl;
    //    } else if (vi1 == 8) {
    //      std::cout << "fi " << fi << " " << j_entry << std::endl;
    //    } else if (vi2 == 8) {
    //      std::cout << "fi " << fi << " " << k_entry << std::endl;
    //    }

    m_divergence(m.getIndex(v0)) +=
        0.5 * ((X_fi * e01) * cotan2 + (X_fi * (-e20)) * cotan1);
    m_divergence(m.getIndex(v1)) +=
        0.5 * ((X_fi * e12) * cotan0 + (X_fi * (-e01)) * cotan2);
    m_divergence(m.getIndex(v2)) +=
        0.5 * ((X_fi * e20) * cotan1 + (X_fi * (-e12)) * cotan0);
  }
  std::cout << "vertex " << debugVi << " " << m_divergence(debugVi)
            << std::endl;

  //  for (size_t vi = 0; vi < m.VN(); vi++) {
  //    outfile << "vertex " << vi << " " << m_divergence(vi) << std::endl;
  //  }
  //  outfile.close();
}

Eigen::VectorXd
GeodesicDistance::solvePhi(const std::unordered_set<VertexIndex> &sourcesVi) {
  Eigen::VectorXd phi = la_cotan.solve(m_divergence);

  Eigen::VectorXd distances(m.VN());
  const double sourceDistance = phi(*sourcesVi.begin(), 0);

  return phi.unaryExpr([&](double x) { return std::abs(-sourceDistance + x); });
  //  for (size_t vi = 0; vi < m.VN(); vi++) {
  //  }

  //    for (size_t vi = 0; vi < m.VN(); vi++) {
  //    std::cout << "Proccessed vertex " << vi << std::endl;
  //    double min_val = (std::numeric_limits<double>::max)();
  //    // go through the distances to the sources and leave the minimum
  //    distance;
  //        for (VertexIndex sourceVi : sourcesVi) {
  //      double new_d = std::abs(-phi(sourceVi, 0) + phi(vi, 0));
  //    if (phi(sourceVi, 0) == phi(vi, 0)) {
  //      min_val = 0.;
  //    }
  //    if (new_d < min_val) {
  //      min_val = new_d;
  //    }
  //        }
  //      distances(vi) = min_val;
  //  }

  //  return distances;
}

Eigen::VectorXd GeodesicDistance::computeGeodesicDistances(
    const std::unordered_set<VertexIndex> &sourcesVi,
    std::unordered_map<VertexIndex, double> &distances) {
  assert(isInitialized);
  //  assert(sourcesVi < m.VN() - 1);

  updateKroneckerDelta(sourcesVi);
  solveCotanLaplace();
  computeUnitGradient();
  computeDivergence();
  return solvePhi(sourcesVi);
}

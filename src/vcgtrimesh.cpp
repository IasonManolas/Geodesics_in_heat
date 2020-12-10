#include "vcgtrimesh.hpp"
#include "wrap/io_trimesh/import_obj.h"
#include "wrap/io_trimesh/import_off.h"
#include <filesystem>

void VCGTriMesh::loadFromPlyFile(const std::string &filename) {
  unsigned int mask = 0;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_VERTCOORD;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_VERTNORMAL;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_VERTCOLOR;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_EDGEINDEX;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_FACEINDEX;
  if (nanoply::NanoPlyWrapper<VCGTriMesh>::LoadModel(
          std::filesystem::absolute(filename).c_str(), *this, mask) != 0) {
    std::cout << "Could not load tri mesh" << std::endl;
  }
  vcg::tri::UpdateTopology<VCGTriMesh>::FaceFace(*this);
  vcg::tri::UpdateTopology<VCGTriMesh>::VertexFace(*this);
  vcg::tri::UpdateNormal<VCGTriMesh>::PerVertexNormalized(*this);
}

Eigen::MatrixX3d VCGTriMesh::getVertices() const {
  Eigen::MatrixX3d vertices(VN(), 3);
  for (size_t vi = 0; vi < VN(); vi++) {
    VCGTriMesh::CoordType vertexCoordinates = vert[vi].cP();
    vertices.row(vi) = vertexCoordinates.ToEigenVector<Eigen::Vector3d>();
  }
  return vertices;
}

Eigen::MatrixX3i VCGTriMesh::getFaces() const {
  Eigen::MatrixX3i faces(FN(), 3);
  for (int fi = 0; fi < FN(); fi++) {
    const VCGTriMesh::FaceType &face = this->face[fi];
    const size_t v0 = vcg::tri::Index<VCGTriMesh>(*this, face.cV(0));
    const size_t v1 = vcg::tri::Index<VCGTriMesh>(*this, face.cV(1));
    const size_t v2 = vcg::tri::Index<VCGTriMesh>(*this, face.cV(2));
    faces.row(fi) = Eigen::Vector3i(v0, v1, v2);
  }
  return faces;
}

bool VCGTriMesh::savePly(const std::string plyFilename) {
  // Load the ply file
  unsigned int mask = 0;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_VERTCOORD;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_VERTCOLOR;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_FACEINDEX;
  mask |= nanoply::NanoPlyWrapper<VCGTriMesh>::IO_FACENORMAL;
  if (nanoply::NanoPlyWrapper<VCGTriMesh>::SaveModel(plyFilename.c_str(), *this,
                                                     mask, false) != 0) {
    return false;
  }
  return true;
}

VCGTriMesh::VCGTriMesh() {}

VCGTriMesh::VCGTriMesh(const std::string &filename) {
  const std::string extension = std::filesystem::path(filename).extension();
  if (extension == ".ply") {
    loadFromPlyFile(filename);
  } else if (extension == ".obj") {
    vcg::tri::io::ImporterOBJ<VCGTriMesh>::Info info;
    vcg::tri::io::ImporterOBJ<VCGTriMesh>::Open(*this, filename.c_str(), info);
  } else if (extension == ".off") {
    vcg::tri::io::ImporterOFF<VCGTriMesh>::Open(*this, filename.c_str());

  } else {
    std::cerr << "Uknown file extension " << extension << ". Could not open "
              << filename << std::endl;
    assert(false);
  }
  vcg::tri::UpdateTopology<VCGTriMesh>::AllocateEdge(*this);
  vcg::tri::UpdateTopology<VCGTriMesh>::FaceFace(*this);
  vcg::tri::UpdateTopology<VCGTriMesh>::VertexFace(*this);
  vcg::tri::UpdateNormal<VCGTriMesh>::PerVertexNormalized(*this);
}

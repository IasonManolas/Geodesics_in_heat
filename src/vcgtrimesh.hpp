#ifndef VCGTRIMESH_HPP
#define VCGTRIMESH_HPP
#include <vcg/complex/complex.h>
#include <wrap/nanoply/include/nanoplyWrapper.hpp>

using VertexIndex = size_t;
class VCGTriMeshVertex;
class VCGTriMeshEdge;
class VCGTriMeshFace;
struct VCGTriMeshUsedTypes
    : public vcg::UsedTypes<vcg::Use<VCGTriMeshVertex>::AsVertexType,
                            vcg::Use<VCGTriMeshEdge>::AsEdgeType,
                            vcg::Use<VCGTriMeshFace>::AsFaceType> {};
class VCGTriMeshVertex
    : public vcg::Vertex<VCGTriMeshUsedTypes, vcg::vertex::Coord3d,
                         vcg::vertex::Normal3d, vcg::vertex::BitFlags,
                         vcg::vertex::Color4b, vcg::vertex::VFAdj> {};
class VCGTriMeshFace
    : public vcg::Face<VCGTriMeshUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj,
                       vcg::face::VertexRef, vcg::face::BitFlags,
                       vcg::face::Normal3d> {};
class VCGTriMeshEdge
    : public vcg::Edge<VCGTriMeshUsedTypes, vcg::edge::VertexRef> {};

class VCGTriMesh : public vcg::tri::TriMesh<std::vector<VCGTriMeshVertex>,
                                            std::vector<VCGTriMeshFace>,
                                            std::vector<VCGTriMeshEdge>> {
public:
  VCGTriMesh();
  VCGTriMesh(const std::string &filename);
  void loadFromPlyFile(const std::string &filename);
  Eigen::MatrixX3d getVertices() const;
  Eigen::MatrixX3i getFaces() const;
  bool savePly(const std::string plyFilename);
  template <typename MeshElement> size_t getIndex(const MeshElement &element) {
    return vcg::tri::Index<VCGTriMesh>(*this, element);
  }
};

#endif // VCGTRIMESH_HPP

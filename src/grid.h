#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cassert>
#include <ranges>

struct v2{
	float x, y;
	v2 operator+(v2 const& r)const {
		return {x+r.x,y+r.y};
	}
	v2 operator-(v2 const& r)const {
		return {x-r.x,y-r.y};
	}
	v2 operator*(float r)const {
		return {x*r,y*r};
	}
	friend std::ostream& operator<<(std::ostream& os,v2 const& v) {
		return os<<v.x<<" "<<v.y;
	}
	bool operator==(v2 const& r)const = default;
};


template<int N> struct always{static constexpr int max(){return N;}};
template<int N> struct upto{static constexpr int max(){return N;}};

using vert_id = int;
using face_id = int;

template<class NFaces> struct vertex;
template<int N> struct vertex<always<N>>{
	v2 pos;
	face_id faces[N];
	vert_id neighbours[N];
	constexpr int getNumFaces() const {return N;}
	void setNumFaces(int n) { assert(n==N); }
};
template<int N> struct vertex<upto<N>> : public vertex<always<N>>{
	unsigned short nFaces = 0;
	int getNumFaces() const {return nFaces;}
	void setNumFaces(int n) { nFaces = n; }
};


template<typename NVertices> struct face;
template<int N> struct face<always<N>>{
	v2 centre;
	vert_id corners[N];
	face_id neighbours[N];
	constexpr int getNumCorners() const {return N;}
	void setNumCorners(int n) { assert(n==N); }
};
template<int N> struct face<upto<N>> : public face<always<N>>{
	unsigned short nCorners = 0;
	int getNumCorners() const {return nCorners;}
	void setNumCorners(int n) { nCorners = n; }
};


template<typename NVertexFaces, typename NFaceVertices>
struct grid {
	using dual_t = grid<NFaceVertices, NVertexFaces>;
	using vertex_t = vertex<NVertexFaces>;
	using face_t = face<NFaceVertices>;
//	struct metaInformation{
//		int complFaceSides[NFaceVertices::max()];
//		int complVertexSides[NVertexFaces::max()];
//	};

	std::vector<vertex_t> vertices;
	std::vector<face_t> faces;
	//metaInformation meta;

	void clear(){
		vertices.clear();
		faces.clear();
	}

	// Ideally, when a face f0 has a face f1 as its d-th neighbour, then f1 has f0 as its phi(d)-th neighbour,
	// where phi is an involutive permutation of {0,...,NFaceVertices-1} and the labels d are circularly arranged (CW or CCW, respectively).
	// Unfortunately, for general grids (especially those with different n-gons) the edges can not be labelled such that a phi exists. int nTotalHoursWastedHere = 16;
	// In these cases, we need to compute the complement index for every edge of every face.
	// It turns out that this is not a performance issue at all. TODO still precompute
	int getComplFace(face_id f, int d) const{
		face_t const& ff = faces[faces[f].neighbours[d]];
		for(int k = 0; k < ff.getNumCorners(); ++k)
			if(ff.neighbours[k] == f)
				return k;
		assert(false);
	}

	int getComplVertex(vert_id v, int d) const{
		vertex_t const& vv = vertices[vertices[v].neighbours[d]];
		for(int k = 0; k < vv.getNumFaces(); ++k)
			if(vv.neighbours[k] == v)
				return k;
		assert(false);
	}

	template<typename CB>
	int forAllFaces(CB&& cb) const {
		int counter{0};
		for (face_t const& f : faces) {
			cb(f);
			counter++;
		}
		return counter;
	}

	void print(std::ostream& os) const {
		for (const face_t& f : faces) {
			for (int k = 0; k < f.getNumCorners(); ++k)
				if(f.corners[k] != -1)
					os << vertices[f.corners[k]].pos << " ";
			if(f.getNumCorners()>0) os << std::endl;
		}
	}
};


// Computes the dual grid. Is an involution.
template<class NVertexFaces, class NFaceVertices>
void computeDualGrid(grid<NVertexFaces,NFaceVertices> const& src, grid<NFaceVertices,NVertexFaces>& dst) {
	dst.clear();
	dst.faces.resize(src.vertices.size());
	dst.vertices.resize(src.faces.size());
	// The old vertices become the new faces
	for (int i = 0; i < src.vertices.size(); ++i) {
		dst.faces[i].centre = src.vertices[i].pos;
		dst.faces[i].setNumCorners(src.vertices[i].getNumFaces());
		for (int j = 0; j < src.vertices[i].getNumFaces(); ++j) {
			dst.faces[i].corners[j] = src.vertices[i].faces[j];
			dst.faces[i].neighbours[j] = src.vertices[i].neighbours[j];
		}
	}
	// The old faces become the new vertices
	for (int i = 0; i < src.faces.size(); ++i) {
		dst.vertices[i].pos = src.faces[i].centre;
		dst.vertices[i].setNumFaces(src.faces[i].getNumCorners());
		for (int j = 0; j < src.faces[i].getNumCorners(); ++j) {
			dst.vertices[i].faces[j] = src.faces[i].corners[j];
			dst.vertices[i].neighbours[j] = src.faces[i].neighbours[j];
		}
	}

//	memcpy(dst.meta.complVertexSides, src.meta.complFaceSides, sizeof(int)*NFaceVertices::max());
//	memcpy(dst.meta.complFaceSides, src.meta.complVertexSides, sizeof(int)*NVertexFaces::max());
}

template<class NVertexFaces, class NFaceVertices>
int checkIntegrity(grid<NVertexFaces,NFaceVertices> const& grid) {
	AutoTimer at(g_timer, _FUNC_);
	int nDefects = 0;
	// check if the face complements are as advertised
	for (int i = 0; i < grid.faces.size(); ++i) {
		const auto& fi = grid.faces[i];
		for (int k = 0; k < fi.getNumCorners(); ++k) {
			auto j = fi.neighbours[k];
			if (j < 0) continue;
			const auto& fj = grid.faces[j];
			if (i == j) ++nDefects;
			if (i != fj.neighbours[grid.getComplFace(i, k)]) ++nDefects;
		}
	}
	std::cout << nDefects << std::endl;
	nDefects = 0;

	// check if the vertex complements are as advertised
	for (int i = 0; i < grid.vertices.size(); ++i) {
		const auto& vi = grid.vertices[i];
		for (int k = 0; k < vi.getNumFaces(); ++k) {
			auto j = vi.neighbours[k];
			if (j < 0) continue;
			const auto& vj = grid.vertices[j];
			if (i == j) ++nDefects;
			if (i != vj.neighbours[grid.getComplVertex(i, k)]) ++nDefects;
		}
	}
	std::cout << nDefects << std::endl;
	return nDefects;
}


template<class Grid>
class GridGenerator {
public:
	using grid_t = Grid;
	virtual void convert(grid_t& g) const = 0;
};


class HexQuadGenerator : public GridGenerator<grid<upto<6>,always<4>>> {

	struct c_vertex : public grid_t::vertex_t {
		v2 force;
		bool fixedX = false;
		bool fixedY = false;
		int nNeighbours = 0;
	};
	struct c_face : public grid_t::face_t {
		bool exists = false;
	};

	std::vector<c_vertex> vertices;
	std::vector<c_face> faces;

	const int mergeProb = 100;
	float totalArea;

	void sortFacesCircularly() {
		for (auto& v : vertices) {
			if(v.getNumFaces() > 2)
			std::sort(v.faces, v.faces + v.getNumFaces(), [&v, this](face_id a, face_id b) {
				auto pa = faces[a].centre - v.pos;
				auto pb = faces[b].centre - v.pos;
				return std::atan2(pa.x,pa.y) < std::atan2(pb.x,pb.y);
			});
		}
	}

public:
	/*  Example for argument n=1:
	 *
	 *      i=0       i=1       i=2
	 *  j=0  0---------2---------4
	 *        \       . \         \
	 *         \    .    \         \
	 *          \ .       \         \
	 *      j=1 10---------12--------14
	 *            \         \       . \
	 *             \         \    .    \
	 *              \         \ .       \
	 *         j=2  20--------22---------24
	 */
	void generate(int n) {
		AutoTimer at(g_timer,_FUNC_);
		n *= 2; // regular hexagons need even n
		faces.reserve(n * 2 * n * 3);
		faces.resize(n * 2 * n);
		int nblocks = faces.size();
		vertices.reserve((n + 1 + n) * (n + 1 + n) + nblocks);// reserve enough for the worst case. we cannot reallocate!
		vertices.resize((n + 1 + n) * (n + 1 + n));

		const float h = std::sqrt(3.f) / 2.f;
		auto getNode = [this, n](int i, int j, int ii = 0, int jj = 0)->c_vertex& {return vertices[(j * 2 + jj) * (n + 1 + n) + (i * 2 + ii)]; };
		auto getNodeIdx = [n](int i, int j, int ii = 0, int jj = 0)->vert_id {return (j * 2 + jj) * (n + 1 + n) + (i * 2 + ii); };
		auto getIdx = [this](c_vertex* v){return std::distance(vertices.data(), v);};
		auto middleNode = [](vert_id v0, vert_id v1)->vert_id {return std::midpoint(v0,v1); }; // the geographic midpoint is also the memory midpoint
		auto isBorder = [n](int i, int j, int ii = 0, int jj = 0) {return i == 0 && ii == 0 || j == 0 && jj == 0 || i == n && ii == 0 || j == n && jj == 0
		                                                           || (i + j) * 2 + ii + jj == n || (i + j) * 2 + ii + jj == n * 3; };
		for (int j = 0; j <= n; ++j) {
			for (int i = 0; i <= n; ++i) {
				auto& nod = getNode(i, j);
				nod.pos = { i + j * 0.5f,j * h };
				nod.fixedX = nod.fixedY = isBorder(i, j);
				if (i > 0) {
					auto& nod = getNode(i, j, -1, 0);
					nod.pos = (getNode(i - 1, j).pos + getNode(i, j).pos) * 0.5f;
					//if (isBorder(i, j)) nod.fixedX = nod.fixedY = true;
					nod.fixedX = nod.fixedY = isBorder(i, j, -1, 0);
				}
				if (j > 0) {
					auto& nod = getNode(i, j, 0, -1);
					nod.pos = (getNode(i, j - 1).pos + getNode(i, j).pos) * 0.5f;
					// if (isBorder(i, j)) nod.fixedX = nod.fixedY = true;
					nod.fixedX = nod.fixedY = isBorder(i, j, 0, -1);
				}
				if (j == 0 || i == 0) continue;

				auto& nod2 = getNode(i, j, -1, -1);
				nod2.pos = (getNode(i, j - 1, -1, 0).pos + getNode(i, j, -1, 0).pos) * 0.5f;
				nod2.fixedX = nod2.fixedY = isBorder(i, j, -1, -1);

				// Create two new districts
				c_face& d0 = faces[((j - 1) * n + (i - 1)) * 2];
				c_face& d1 = faces[((j - 1) * n + (i - 1)) * 2 + 1];

				if (i + j > n / 2 + 1 && i + j < n + n / 2 + 2) { // master
					d0.exists = true;
					d0.corners[0] = getNodeIdx(i - 1, j - 1);
					d0.corners[1] = getNodeIdx(i, j - 1);
					d0.corners[3] = getNodeIdx(i - 1, j);
					d0.corners[2] = d0.corners[1];
				}
				if (i + j > n / 2 && i + j < n + n / 2 + 1) { // slave
					d1.exists = true;
					d1.corners[0] = getNodeIdx(i, j - 1);
					d1.corners[2] = getNodeIdx(i, j);
					d1.corners[3] = getNodeIdx(i - 1, j);
					d1.corners[1] = d1.corners[0];
				}

				// Sometimes the left district merges with its outer neighbours, or the two new districts merge together:		
				if (d0.exists && rand() % 100 < mergeProb) {
					switch (rand() % 3) {
					case 0:
						if (j >= 2) { // merge up
							c_face& da = faces[((j - 2) * n + (i - 1)) * 2 + 1];
							if (da.exists) {
								d0.corners[1] = da.corners[1];
								da.exists = false;
								break;
							}
						}
					case 1:
						if (i >= 2) { // merge left
							c_face& da = faces[((j - 1) * n + (i - 2)) * 2 + 1];
							if (da.exists) {
								d0.corners[2] = d0.corners[3];
								d0.corners[3] = da.corners[3];
								da.exists = false;
							}
						}
						break;
					case 2:
						if (d0.exists && d1.exists) { // merge self
							d0.corners[2] = d1.corners[2];
							d1.exists = false;
						}
					}
				}
				// The left district is definitely finished now -> create 3 or 4 subcells:
				auto createSubcells = [&, this](c_face& d) {
					if (!d.exists) return;
					v2 center{ 0, 0 }; bool tri = false;
					for (int k = 0; k < 4; ++k) {
						const auto& p = vertices[d.corners[k]].pos;
						const auto& p0 = vertices[d.corners[(k + 3) % 4]].pos;
						if (p == p0)
							tri = true;
						else
							center = center + p;
					}
					center = center * (1.f / (tri ? 3 : 4));
					if (tri) {
						c_vertex& centerNode = vertices.emplace_back();
						centerNode.pos = center;
						const int idx[3] = { 0,2,3 };
						for (int k = 0; k < 3; ++k) {
							c_face& cn = faces.emplace_back();
							face_id cnIdx = faces.size()-1;
							cn.exists = true;
							cn.corners[0] = getIdx(&centerNode);
							cn.corners[1] = middleNode(d.corners[idx[(k + 2) % 3]], d.corners[idx[k]]);
							cn.corners[2] = d.corners[idx[k]];
							cn.corners[3] = middleNode(d.corners[idx[k]], d.corners[idx[(k + 1) % 3]]);
							for (int kk = 0; kk < 4; ++kk)
								vertices[cn.corners[kk]].faces[vertices[cn.corners[kk]].nFaces++] = cnIdx;
						}
					}
					else {
						vert_id centerNodeIdx = middleNode(d.corners[0], d.corners[2]);
						c_vertex& centerNode = vertices[centerNodeIdx];
						for (int k = 0; k < 4; ++k) {
							c_face& cn = faces.emplace_back();
							face_id cnIdx = faces.size()-1;
							cn.exists = true;
							cn.corners[0] = centerNodeIdx;
							cn.corners[1] = middleNode(d.corners[(k + 3) % 4], d.corners[k]);
							cn.corners[2] = d.corners[k];
							cn.corners[3] = middleNode(d.corners[k], d.corners[(k + 1) % 4]);
							for (int kk = 0; kk < 4; ++kk)
								vertices[cn.corners[kk]].faces[vertices[cn.corners[kk]].nFaces++] = cnIdx;
						}
					}
					d.exists = false;
				};

				createSubcells(d0); // definitely finished
				if (i > 1 && j > 1)
					createSubcells(faces[((j - 2) * n + (i - 2)) * 2 + 1]); // completely enclosed now
				if (i > 1 && j == n)
					createSubcells(faces[((j - 1) * n + (i - 2)) * 2 + 1]); // enclosed at bottom
				if (i == n && j > 1)
					createSubcells(faces[((j - 2) * n + (i - 1)) * 2 + 1]); // enclosed at right side
			}
		}

		for(auto& v: vertices) {
			for (int k = 0; k < 6; ++k) {
				v.neighbours[k] = -1;
			}
		}

		// Use corner information to find the neighbours of each face:
		for (c_face& f0 : faces) {
			if (!f0.exists) continue;
			for (int d0 = 0; d0 < 4; ++d0) {
				f0.neighbours[d0] = -1;
				auto& v = vertices[f0.corners[d0]];
				for (int k = 0; k < v.nFaces; ++k) {
					face_id f1 = v.faces[k];
					if (&faces[f1] == &f0) continue;
					for (int d1 = 0; d1 < 4; ++d1) {
						vert_id v0 = faces[f1].corners[d1], v1 = f0.corners[(d0 + 1) % 4];
						if (v0 == v1) {
							f0.neighbours[d0] = f1;
							v.neighbours[v.nNeighbours++] = v1;
							goto next_side;
						}
					}
				}
				next_side:
				continue;
			}
		}

		std::map<int,int> usedVerts;
		for (c_face& fa : faces) {
			if(fa.exists)
				for(int i = 0; i< fa.getNumCorners();++i)
					usedVerts[fa.corners[i]]++;
		}
		std::map<int,int> counter;
		for(auto& v: vertices) {
			int i =&v-vertices.data() ;
			if(usedVerts[i] > 0)
				counter[v.nFaces]++;
		}

		// Compute face centres
		for (c_face& f : faces) {
			f.centre = {0, 0};
			int nCorners = 0;
			for (int i = 0; i < f.getNumCorners(); ++i)
				if (f.corners[i] >= 0) {
					f.centre = f.centre + vertices[f.corners[i]].pos;
					nCorners++;
				}
			f.centre = f.centre * (1.f / nCorners);
		}

		// Make sure that the faces of each vertex are ordered circularly. This will make the dual grid trivially drawable.
		sortFacesCircularly();

		totalArea = std::sqrt(3.f) * 3 / 2 * n * n / 4;
	}
	void relax(int nRelaxIterations){

		AutoTimer at(g_timer,_FUNC_);
		int ncellstot = 0;
		for (c_face& c : faces) if(c.exists) ncellstot++;
		float avgArea = totalArea / ncellstot;
		for (int iter = 0; iter < nRelaxIterations; iter++) {
			for (c_face& c : faces) {
				if (!c.exists) continue;
				v2 center = { 0,0 };
				float area = 0;
				for (int i = 0; i < 4; ++i) {
					const auto& p = vertices[c.corners[i]].pos;
					const auto& p0 = vertices[c.corners[(i + 3) % 4]].pos;
					center = center + p;
					area += (p0.x * p.y - p0.y * p.x) / 2;
				}
				center.x /= 4;
				center.y /= 4;
				area = (abs(area) - avgArea) / avgArea;
				v2 force = { 0,0 };
				for (int i = 0; i < 4; ++i) {
					force = force + (vertices[c.corners[i]].pos - center);
					force = { -force.y, force.x };
				}
				for (int i = 0; i < 4; ++i) {
					vertices[c.corners[i]].force = vertices[c.corners[i]].force + (center - vertices[c.corners[i]].pos) * (1 + area) + force * 0.25;
					force = { -force.y, force.x };
				}
			}
			for (int j = 0; j < vertices.size(); ++j) {
				c_vertex& k = vertices[j];
				if (!k.fixedX)
					k.pos.x = k.pos.x + k.force.x * 0.2f;
				if (!k.fixedY)
					k.pos.y = k.pos.y + k.force.y * 0.2f;
				k.force = { 0,0 };
			}
		}
	}


	void convert(grid_t& g) const override {
		// We want to filter the non-existent faces.
		// Create a map from old to new face indices:
		std::vector<face_id> ids(faces.size(), 1);
		for(int i = 0; i< faces.size(); ++i)
			ids[i] = faces[i].exists;
		std::exclusive_scan(ids.begin(), ids.end(), ids.begin(), 0);

		g.clear();
		g.vertices.resize(vertices.size());
		g.faces.reserve(ids.back());

		auto convertVertex = [&ids](c_vertex const& cv){
			grid_t::vertex_t v = cv;
			for(int i =0; i< v.nFaces; ++i)
				v.faces[i] = ids[v.faces[i]];
			return v;
		};
		auto convertFace = [&ids, this](c_face const& cf) {
			face f = cf;
			for (int i = 0; i < cf.getNumCorners(); ++i){
				if (f.neighbours[i] >= 0) {
					f.neighbours[i] = ids[f.neighbours[i]];
				}
			}
			return f;
		};

		std::transform(vertices.begin(), vertices.end(), g.vertices.begin(), convertVertex);

		for(auto a : faces | std::views::filter([](c_face const& cf){return cf.exists;}) | std::views::transform(convertFace))
			g.faces.push_back(a);


//		g.meta.complFaceSides[0] = 3;
//		g.meta.complFaceSides[1] = 2;
//		g.meta.complFaceSides[2] = 1;
//		g.meta.complFaceSides[3] = 0;
//
//		g.meta.complVertexSides[0] = 0;
//		g.meta.complVertexSides[1] = 0;
//		g.meta.complVertexSides[2] = 0;
//		g.meta.complVertexSides[3] = 0;
//		g.meta.complVertexSides[4] = 0;
//		g.meta.complVertexSides[5] = 0;
	}

	void print() const{
		std::ofstream file("debug.txt");
		for(int i=0; auto const& v :vertices){
			file << i << " " << v.pos;
			for(int k = 0; k< v.getNumFaces(); ++k)
				file << " " << v.neighbours[k];
			file << std::endl;
			++i;
		}

		file.close();
	}

};

class SquareGenerator : public GridGenerator<grid<always<4>,always<4>>> {

	struct c_vertex : public grid_t::vertex_t {
	};
	struct c_face : public grid_t::face_t {
	};

	std::vector<c_vertex> vertices;
	std::vector<c_face> faces;

public:
	/*  Example for argument n=1:
	 *
	 *      i=0       i=1       i=2
	 *  j=0  0---------1---------2
	 *       |         |         |
	 *       |         |         |
	 *       |         |         |
	 *  j=1  3---------4---------5
	 *       |         |         |
	 *       |         |         |
	 *       |         |         |
	 *  j=2  6---------7---------8
	 */
	void generate(int n) {
		AutoTimer at(g_timer, _FUNC_);

		faces.resize(n * n);
		vertices.resize((n + 1) * (n + 1));

		for (auto& a: vertices)
			for (int i = 0; i < 4; ++i)
				a.faces[i] = -1;

		auto getNode = [this, n](int i, int j) -> c_vertex& { return vertices[j * (n + 1) + i]; };
		auto getNodeIdx = [n](int i, int j) -> vert_id { return j * (n + 1) + i; };
		for (int j = 0; j <= n; ++j) {
			for (int i = 0; i <= n; ++i) {
				auto& nod = getNode(i, j);
				nod.pos = {(float) i, (float) j};

				if (j == 0 || i == 0) continue;

				// Create a new face
				face_id fid = (i - 1) + n * (j - 1);
				c_face& f0 = faces[fid];
				f0.corners[0] = getNodeIdx(i - 1, j - 1);
				f0.corners[1] = getNodeIdx(i, j - 1);
				f0.corners[3] = getNodeIdx(i - 1, j);
				f0.corners[2] = getNodeIdx(i, j);
				for (int kk = 0; kk < 4; ++kk)
					vertices[f0.corners[kk]].faces[(kk + 2) % 4] = fid;
			}
		}


		// Use node information to find the neighbours of each cell:
		for (c_face& f : faces) {
			for (int kk = 0; kk < f.getNumCorners(); ++kk) {
				f.neighbours[kk] = -1;
				const auto& v = vertices[f.corners[kk]];
				for (int xy = 0; xy < v.getNumFaces(); ++xy) {
					if (v.faces[xy] == -1) continue;
					const c_face& otherface = faces[v.faces[xy]];
					if (&otherface == &f) continue;
					bool touches = false;
					for (int kk2 = 0; kk2 < 4; ++kk2) {
						if (otherface.corners[kk2] == f.corners[(kk + 1) % 4]) {
							touches = true;
							break;
						}
					}
					if (touches) {
						f.neighbours[kk] = &otherface - faces.data();
						break;
					}
				}
			}
		}

		// Compute face centres
		for (c_face& f : faces) {
			f.centre = {0, 0};
			int nCorners = 0;
			for (int i = 0; i < f.getNumCorners(); ++i)
				if (f.corners[i] >= 0) {
					f.centre = f.centre + vertices[f.corners[i]].pos;
					nCorners++;
				}
			f.centre = f.centre * (1.f / nCorners);
		}
	}


	void convert(grid_t& g) const override {
		// We want to filter the non-existent faces.
		// Create a map from old to new face indices:

		g.clear();
		g.vertices.resize(vertices.size());
		g.faces.reserve(faces.size());

		auto convertVertex = [](c_vertex const& cv){
			grid_t::vertex_t v = cv;
			for(int i =0; i < v.getNumFaces(); ++i)
				v.faces[i] = v.faces[i];
			return v;
		};
		auto convertFace = [this](c_face const& cf){
			grid_t::face_t f = cf;
			return f;
		};


		std::transform(vertices.begin(), vertices.end(), g.vertices.begin(), convertVertex);
		//	std::transform(faces.begin(), faces.end(), g.faces.begin(), convertFace);

		for(auto a : faces | std::views::transform(convertFace))
			g.faces.push_back(a);

//		g.meta.complFaceSides[0] = 2;
//		g.meta.complFaceSides[1] = 3;
//		g.meta.complFaceSides[2] = 0;
//		g.meta.complFaceSides[3] = 1;
//
//		g.meta.complVertexSides[0] = 0;
//		g.meta.complVertexSides[1] = 0;
//		g.meta.complVertexSides[2] = 0;
//		g.meta.complVertexSides[3] = 0;
	}

};

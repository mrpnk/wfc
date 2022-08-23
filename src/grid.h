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


//
//template<int N>
//struct always;
//
//template<int N>
//struct upTo;
//
//
//template<typename T>
//struct vertex;
//
//
//template<int N>
//struct vertex<always<N>>{
//	v2 pos;
//	int neighbours[N];
//};
//
//template<int N>
//struct vertex<upTo<N>> : private vertex<always<N>>{
//	int num;
//};
//
//vertex<always<4>> sdf;
//vertex<upTo<4>> sddf;
//

using vert_id = int;
using face_id = int;

struct vertex {
	v2 pos;
	unsigned short nAdjacentFaces = 0;
	face_id adjacentFaces[6];
};
struct face {
	vert_id corners[4];
	face_id nn[4];
	mutable void* data = nullptr;
	mutable bool touched = false;
};

struct Grid {
	std::vector<vertex> vertices;
	std::vector<face> faces;


	template<typename CB>
	int forAllCells(CB&& cb) const {
		int counter = 0;
		for (face const& c : faces) {
			//if (!c.exists) continue;
			cb(c);
			counter++;
		}
		return counter;
	}

	void print(std::ostream& os) {
		for (face& c : faces) {
			//if (!c.exists) continue;
			for (int k = 0; k < 4; ++k)
				os << vertices[c.corners[k]].pos << " ";
			os << 0 << " ";
			os << std::endl;
		}
	}
};

class GridGenerator{
	struct c_vertex : public vertex{
		v2 force;
		bool fixedX = false;
		bool fixedY = false;
	};
	struct c_face : public face {
		bool exists = false;
	};

	std::vector<c_vertex> vertices;
	std::vector<c_face> faces;

	const int mergeProb = 100;
	float totalArea;


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
	void construct(int n) {
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
								vertices[cn.corners[kk]].adjacentFaces[vertices[cn.corners[kk]].nAdjacentFaces++] = cnIdx;
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
								vertices[cn.corners[kk]].adjacentFaces[vertices[cn.corners[kk]].nAdjacentFaces++] = cnIdx;
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


		// Use node information to find the neighbours of each cell:
		int haven[5] = { 0,0,0,0,0 };
		for (c_face& fa : faces) {
			if (!fa.exists) continue;
			int nnn = 0;
			for (int kk = 0; kk < 4; ++kk) {
				fa.nn[kk] = -1;
				const auto& v = vertices[fa.corners[kk]];
				for (int xy = 0; xy < v.nAdjacentFaces; ++xy) {
					face_id otherface = v.adjacentFaces[xy];
					if (&faces[otherface] == &fa) continue;
					bool touches = false;
					for (int kk2 = 0; kk2 < 4; ++kk2) {
						if (faces[otherface].corners[kk2] == fa.corners[(kk + 1) % 4]) {
							touches = true;
							break;
						}
					}
					if (touches){
						nnn++;
						fa.nn[kk] = otherface;
						break;
					}
				}
			}
			haven[nnn]++;
		}
//		for (int i = 0; i <= 4; ++i)
//			std::cout << "have " << i << ": " << haven[i] << std::endl;

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


	void convert(Grid& g) const {
		// We want to filter the non-existent faces.
		// Create a map from old to new face indices:
		std::vector<face_id> ids(faces.size(), 1);
		for(int i = 0; i< faces.size(); ++i)
			ids[i] = faces[i].exists;
		std::exclusive_scan(ids.begin(), ids.end(), ids.begin(), 0);

		g.vertices.resize(vertices.size());
		g.faces.reserve(ids.back());

		auto convertVertex = [&ids](c_vertex const& cv){
			vertex v = cv;
			for(int i =0; i< v.nAdjacentFaces; ++i)
				v.adjacentFaces[i] = ids[v.adjacentFaces[i]];
			return v;
		};
		auto convertFace = [&ids](c_face const& cf){
			face f = cf;
			for(int i =0; i < 4; ++i)
				if(f.nn[i] >=0) f.nn[i] = ids[f.nn[i]];
			return f;
		};



		std::transform(vertices.begin(), vertices.end(), g.vertices.begin(), convertVertex);
	//	std::transform(faces.begin(), faces.end(), g.faces.begin(), convertFace);

		for(auto a : faces | std::views::filter([](c_face const& cf){return cf.exists;}) | std::views::transform(convertFace))
			g.faces.push_back(a);

		//c_face test{.exists = true, .pos = {23,345}};

	}

};

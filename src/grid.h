#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cassert>

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

class Grid {
public:
	struct cell;
	struct vertex {
		v2 pos;
		v2 force;
		bool fixedX = false;
		bool fixedY = false;
		int ncontained = 0;
		cell* contained[6] = { 0,0,0,0,0,0 };
	};
	struct cell {
		bool exists = false;
		vertex* corners[4];
		cell* nn[4];
		mutable void* data = nullptr;
		mutable bool touched = false;
	};
private:
	std::vector<vertex> vertices;
	std::vector<cell> cells;
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
	 *      j=1  5---------7---------9
	 *            \         \       . \
	 *             \         \    .    \
	 *              \         \ .       \
	 *         j=2  10--------12---------14
	 *
	 */
	void init(int n) {
		AutoTimer at(g_timer,_FUNC_);
		n *= 2; // regular hexagons need even n
		cells.reserve(n * 2 * n * 3);
		cells.resize(n * 2 * n);
		int nblocks = cells.size();
		vertices.reserve((n + 1 + n) * (n + 1 + n) + nblocks);// reserve enough for the worst case. we cannot reallocate!
		vertices.resize((n + 1 + n) * (n + 1 + n));

		float h = std::sqrt(3.f) / 2.f;
		auto getNode = [this, n](int i, int j, int ii = 0, int jj = 0)->vertex& {return vertices[(j * 2 + jj) * (n + 1 + n) + (i * 2 + ii)]; };
		auto middleNode = [](vertex* nod0, vertex* nod1)->vertex& {return *(nod0 + std::distance(nod0, nod1) / 2); };
		auto isBorder = [n](int i, int j, int ii = 0, int jj = 0) {return i == 0 || j == 0 || i == n || j == n || (i + j) * 2 + ii + jj == n || (i + j) * 2 + ii + jj == n * 3; };
		for (int j = 0; j <= n; ++j) {
			for (int i = 0; i <= n; ++i) {
				vertex& nod = getNode(i, j);
				nod.pos = { i + j * 0.5f,j * h };
				nod.fixedX = nod.fixedY = isBorder(i, j);
				if (i > 0) {
					vertex& nod = getNode(i, j, -1, 0);
					nod.pos = (getNode(i - 1, j).pos + getNode(i, j).pos) * 0.5f;
					if (isBorder(i, j)) nod.fixedX = nod.fixedY = true;
				}
				if (j > 0) {
					vertex& nod = getNode(i, j, 0, -1);
					nod.pos = (getNode(i, j - 1).pos + getNode(i, j).pos) * 0.5f;
					if (isBorder(i, j)) nod.fixedX = nod.fixedY = true;
				}
				if (j == 0 || i == 0) continue;

				vertex& nod2 = getNode(i, j, -1, -1);
				nod2.pos = (getNode(i, j - 1, -1, 0).pos + getNode(i, j, -1, 0).pos) * 0.5f;
				if (isBorder(i, j, -1, -1)) nod2.fixedX = nod2.fixedY = true;

				// Create two new districts
				cell& d0 = cells[((j - 1) * n + (i - 1)) * 2];
				cell& d1 = cells[((j - 1) * n + (i - 1)) * 2 + 1];

				if (i + j > n / 2 + 1 && i + j < n + n / 2 + 2) { // master
					d0.exists = true;
					d0.corners[0] = &getNode(i - 1, j - 1);
					d0.corners[1] = &getNode(i, j - 1);
					d0.corners[3] = &getNode(i - 1, j);
					d0.corners[2] = d0.corners[1];
				}
				if (i + j > n / 2 && i + j < n + n / 2 + 1) { // slave
					d1.exists = true;
					d1.corners[0] = &getNode(i, j - 1);
					d1.corners[2] = &getNode(i, j);
					d1.corners[3] = &getNode(i - 1, j);
					d1.corners[1] = d1.corners[0];
				}

				// Sometimes the left district merges with its outer neighbours, or the two new districts merge together:		
				if (d0.exists && rand() % 100 < mergeProb) {
					switch (rand() % 3) {
					case 0:
						if (j >= 2) { // merge up
							cell& da = cells[((j - 2) * n + (i - 1)) * 2 + 1];
							if (da.exists) {
								d0.corners[1] = da.corners[1];
								da.exists = false;
								break;
							}
						}
					case 1:
						if (i >= 2) { // merge left
							cell& da = cells[((j - 1) * n + (i - 2)) * 2 + 1];
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
				auto createSubcells = [&, this](cell& d) {
					if (!d.exists) return;
					v2 center = { 0,0 }; bool tri = false;
					for (int k = 0; k < 4; ++k) {
						auto& p = d.corners[k]->pos;
						auto& p0 = d.corners[(k + 3) % 4]->pos;
						if (p == p0)
							tri = true;
						else
							center = center + p;
					}
					center = center * (1.f / (tri ? 3 : 4));
					if (tri) {
						vertex& centerNode = vertices.emplace_back(vertex{ .pos = center });
						int idx[3] = { 0,2,3 };
						for (int k = 0; k < 3; ++k) {
							cell& cn = cells.emplace_back();
							cn.exists = true;
							cn.corners[0] = &centerNode;
							cn.corners[1] = &middleNode(d.corners[idx[(k + 2) % 3]], d.corners[idx[k]]);
							cn.corners[2] = d.corners[idx[k]];
							cn.corners[3] = &middleNode(d.corners[idx[k]], d.corners[idx[(k + 1) % 3]]);
							for (int kk = 0; kk < 4; ++kk)
								cn.corners[kk]->contained[cn.corners[kk]->ncontained++] = &cn;
						}
					}
					else {
						vertex* centerNode = &middleNode(d.corners[0], d.corners[2]);
						for (int k = 0; k < 4; ++k) {
							cell& cn = cells.emplace_back();
							cn.exists = true;
							cn.corners[0] = centerNode;
							cn.corners[1] = &middleNode(d.corners[(k + 3) % 4], d.corners[k]);
							cn.corners[2] = d.corners[k];
							cn.corners[3] = &middleNode(d.corners[k], d.corners[(k + 1) % 4]);
							for (int kk = 0; kk < 4; ++kk)
								cn.corners[kk]->contained[cn.corners[kk]->ncontained++] = &cn;
						}
					}
					d.exists = false;
				};

				createSubcells(d0); // definitely finished
				if (i > 1 && j > 1)
					createSubcells(cells[((j - 2) * n + (i - 2)) * 2 + 1]); // completely enclosed now
				if (i > 1 && j == n)
					createSubcells(cells[((j - 1) * n + (i - 2)) * 2 + 1]); // enclosed at bottom
				if (i == n && j > 1)
					createSubcells(cells[((j - 2) * n + (i - 1)) * 2 + 1]); // enclosed at right side
			}
		}


		// Use node information to find the neighbours of each cell:
		int haven[5] = { 0,0,0,0,0 };
		for (cell& face : cells) {
			if (!face.exists) continue;
			int nnn = 0;
			for (int kk = 0; kk < 4; ++kk) {
				face.nn[kk] = nullptr;
				for (int xy = 0; xy < face.corners[kk]->ncontained; ++xy) {
					cell* otherface = face.corners[kk]->contained[xy];
					if (otherface == &face) continue;
					bool touches = false;
					for (int kk2 = 0; kk2 < 4; ++kk2) {
						if (otherface->corners[kk2] == face.corners[(kk + 1) % 4]) {
							touches = true;
							break;
						}
					}
					if (touches)
					{
						nnn++;
						face.nn[kk] = otherface;
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
		for (cell& c : cells)if (c.exists)ncellstot++;
		float avgArea = totalArea / ncellstot;
		for (int iter = 0; iter < nRelaxIterations; iter++) {
			for (cell& c : cells) {
				if (!c.exists) continue;
				v2 center = { 0,0 };
				float area = 0;
				for (int i = 0; i < 4; ++i) {
					auto p = c.corners[i]->pos;
					auto p0 = c.corners[(i + 3) % 4]->pos;
					center = center + p;
					area += (p0.x * p.y - p0.y * p.x) / 2;
				}
				center.x /= 4;
				center.y /= 4;
				area = (abs(area) - avgArea) / avgArea;
				v2 force = { 0,0 };
				for (int i = 0; i < 4; ++i) {
					force = force + (c.corners[i]->pos - center);
					force = { -force.y, force.x };
				}
				for (int i = 0; i < 4; ++i) {
					c.corners[i]->force = c.corners[i]->force + (center - c.corners[i]->pos) * (1 + area) + force * 0.25;
					force = { -force.y, force.x };
				}
			}
			for (int j = 0; j < vertices.size(); ++j) {
				vertex& k = vertices[j];
				if (!k.fixedX)
					k.pos.x = k.pos.x + k.force.x * 0.2f;
				if (!k.fixedY)
					k.pos.y = k.pos.y + k.force.y * 0.2f;
				k.force = { 0,0 };
			}
		}
	}

	void print(std::ostream& os) {
		for (cell& c : cells) {
			if (!c.exists) continue;
			for (int k = 0; k < 4; ++k)
				os << c.corners[k]->pos << " ";
			os << 0 << " ";
			os << std::endl;
		}
	}

	template<typename CB>
	int forAllCells(CB&& cb) const {
		int counter = 0;
		for (cell const& c : cells) {
			if (!c.exists) continue;
			cb(c);
			counter++;
		}
		return counter;
	}

};

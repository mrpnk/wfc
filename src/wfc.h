#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stack>
#include <queue>
#include <set>
#include <cassert>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <mutex>

#include "grid.h"

namespace wfc
{
	template<typename T>
	class uniqueStack {
		std::stack<T> st;
		std::set<T> se;
		mutable std::mutex mtx;
	public:
		bool push(T i) {
			const std::lock_guard lock(mtx);
			if (se.contains(i))
				return false;
			st.push(i);
			se.insert(i);
			return true;
		}
		bool empty() { return st.empty(); }
		T pop() {
			auto temp = st.top();
			st.pop();
			se.erase(temp);
			return temp;
		}
		int size() const {
			assert(st.size() == se.size());
			return st.size();
		}
	};

	struct colour {
		float r, g, b;
		bool operator<(colour const& c) const { return r < c.r ? true : r > c.r ? false : g < c.g ? true : g > c.g ? false : b < c.b; }
		bool operator==(colour const& c) const = default;
		friend std::ostream& operator<<(std::ostream& os, colour const& col) { return os << col.r << " " << col.g << " " << col.b; }
	};

	template<int m> struct segment {
		colour cols[m * m];
	};


	using superposition = std::vector<int>;
	using wavefunction = std::vector<superposition>;

	static const int m = 6; //tile size per axis in pixels
	static const int k = 5; //source texture size per axis in tiles
	static const bool MIRROR = 1;

	using tt = unsigned short;
	using seg = segment<m>;
	struct option {
		tt modIdx;
		tt rot;
		bool mirrorX;
		bool operator==(option const&)const = default;
	};


	struct gridState{
		std::vector<seg> modus;
		std::vector<option> allOptions;
		std::vector<bool> fitlus;
		wavefunction gridwave;
		std::vector<int> allSP;

		gridState(){
			const int nmodus = k * k;
			const int noptions = nmodus * 4 * (MIRROR + 1);
			modus.resize(nmodus);
			allOptions.resize(noptions);
			fitlus.resize((nmodus * 4) * (nmodus * 4) * (MIRROR + 1));
		}
		void load(){
			std::ifstream infile("modules.txt", std::ios::binary);
			colour imagedata[m * m * k * k];
			infile.read((char*)(imagedata), sizeof(imagedata));
			infile.close();
			for (int i = 0; i < k; ++i) {
				for (int j = 0; j < k; ++j) {
					for (int a = 0; a < m; ++a)
						for (int b = 0; b < m; ++b)
							modus[(i * k + j)].cols[a * m + b] = imagedata[m * m * k * i + m * j + m * k * a + b];
				}
			}
			std::cout << "Extracted modules." << std::endl;
			for (tt i = 0; i < modus.size(); ++i)
				for (tt r = 0; r < 4; ++r)
					for (tt mirr = 0; mirr <= MIRROR; ++mirr)
						allOptions[(i * 4 + r) * (MIRROR + 1) + mirr] = { i,r,(bool)mirr };
			allSP.resize(allOptions.size()); iota(begin(allSP), end(allSP), 0);
			computeFitLU();
			srand(77);//std::chrono::system_clock::now().time_since_epoch().count());
		}
		void initWithGrid(Grid const& g){
			// Create a superposition for each cell in the grid:
			gridwave.resize(g.forAllCells([](auto) {}), allSP);
			g.forAllCells([&, counter = 0](Grid::cell const& c)mutable{ c.data = &gridwave[counter++]; });
		}

		Grid::cell const* getPos(Grid::cell const* p, int d) const {
			return p->nn[d];
		}
		inline int faceElemIdx(int side, int i) const {
			int indices[4] = { i,m - 1 + i * m,m * m - 1 - i,(m - 1 - i) * m };
			return indices[side];
		}

		bool fit(seg const& mod1, int s1, seg const& mod2, int s2, bool mirror) const {
			for (int i = 0; i < m; ++i)
				if (mod1.cols[faceElemIdx(s1, i)] != mod2.cols[faceElemIdx(s2, mirror ? i : m - 1 - i)])
					return false;
			return true;
		}
		void computeFitLU() {
			for (tt mod1 = 0; mod1 < modus.size(); ++mod1) {
				for (tt s1 = 0; s1 < 4; ++s1) {
					for (tt mod2 = 0; mod2 < modus.size(); ++mod2) {
						for (tt s2 = 0; s2 < 4; ++s2) {
							for (tt mirrorX = 0; mirrorX <= MIRROR; ++mirrorX)
								fitlus[(((mod1 * 4 + s1) * modus.size() + mod2) * 4 + s2) * (MIRROR + 1) + mirrorX] = fit(modus[mod1], s1, modus[mod2], s2, (bool)mirrorX);
						}
					}
				}
			}
		}
		inline bool fitlu(tt mod1, tt s1, tt mod2, tt s2, bool mirror) const {
			return fitlus[(((mod1 * 4 + s1) * modus.size() + mod2) * 4 + s2) * (MIRROR + 1) + mirror];
		}
		inline bool fitlu(tt mod1, tt rot1, bool mirror1, tt side1, tt mod2, tt rot2, bool mirror2) const {
			int s1 = (side1 - rot1 + 4) % 4, s2 = (3 - side1 - rot2 + 4) % 4;
			if (mirror1) { s1 = (4 - s1) % 4; }
			if (mirror2) { s2 = (4 - s2) % 4; }
			return fitlu(mod1, s1, mod2, s2, mirror1 ^ mirror2);
		}
		inline bool isPossibleNeighbour(superposition const& sp, tt dir, option const& other_opt) const {
			for (int i : sp)
				if (const option& own_opt = allOptions[i];
						fitlu(own_opt.modIdx, own_opt.rot, own_opt.mirrorX, dir,
						      other_opt.modIdx, other_opt.rot, other_opt.mirrorX))
					return true;
			return false;
		}
		bool isCollapsed() const {
			for (auto& a : gridwave) if (a.size() == 0) return true; else if (a.size() > 1) return false;
			return true;
		}
		bool isStuck() const {
			for (auto& a : gridwave) if (a.size() == 0) return true;
			return false;
		}
		Grid::cell const* getLowestEntropyPos(Grid const& g) const {
			std::vector<Grid::cell const*> minPosis; int minEntropy = std::numeric_limits<int>::max();
			g.forAllCells([&](Grid::cell const& c) {
				superposition& sp = *static_cast<superposition*>(c.data);
				if (sp.size() > 1) {
					if (sp.size() < minEntropy) minPosis = { &c }, minEntropy = sp.size();
					else if (sp.size() == minEntropy) minPosis.push_back(&c);
				}
			});
			return minPosis[rand() % minPosis.size()];
		}
		bool isAnyMirrored() const {
			for (auto& a : gridwave)
				for (auto& b : a)
					if (allOptions[b].mirrorX)
						return true;
			return false;
		}
		void printCanvas(std::ostream& os, std::string sep = "") const {
			for (auto& sp : gridwave) {
				if (sp.size() != 1)
					os << "-1 0 0" << sep;
				else {
					auto opt = allOptions[sp[0]];
					os << opt.modIdx << " " << opt.rot << " " << opt.mirrorX << sep;
				}
				os << std::endl;
			}
		}
	};


	class WaveFunctionCollapser
	{
	public:
		void solve(Grid const& g, gridState& gs) {
			auto t0 = std::chrono::high_resolution_clock::now();
			while (!gs.isCollapsed()) {
				auto p = gs.getLowestEntropyPos(g);
				collapse(gs, p);
				propagate(gs, p);
			}
			std::cout << "WFC done in " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << " ms" << std::endl;
			std::cout << "isAnyMirrored = " << gs.isAnyMirrored() << std::endl;
			std::cout << "isStuck = " << gs.isStuck() << std::endl;
		}

	private:
		seg mirrorX(seg const& mod, bool mirror) const {
			if (!mirror)return mod;
			seg ret;
			for (int a = 0; a < m; ++a) {
				for (int b = 0; b < m; ++b) {
					ret.cols[a * m + b] = mod.cols[a * m + m - 1 - b];
				}
			}
			return ret;
		}
		seg rotate(seg const& mod, int rot) const {
			if (rot == 0) return mod;
			seg ret;
			for (int a = 0; a < m; ++a) {
				for (int b = 0; b < m; ++b) {
					ret.cols[a * m + b] = mod.cols[((m - 1) - b) * m + a];
				}
			}
			return rotate(ret, rot - 1);
		}


		void propagate(gridState& gs, Grid::cell const* pos) {
			uniqueStack<Grid::cell const*> jobs;
			jobs.push(pos);
			while (!jobs.empty()) {
				Grid::cell const* p = jobs.pop();
				superposition& sp = *static_cast<superposition*>(p->data);
				for (int d = 0; d < 4; ++d) {
					if (Grid::cell const* nnp = gs.getPos(p, d); nnp != nullptr) {
						superposition& nnsp = *static_cast<superposition*>(nnp->data);
						if (nnsp.size() <= 1) continue;
						if (std::erase_if(nnsp, [&](int& other_opt) {return !gs.isPossibleNeighbour(sp, d, gs.allOptions[other_opt]); }))
							jobs.push(nnp);
					}
				}
			}
		}
		void collapse(gridState& gs, Grid::cell const* pos) {
			superposition& sp = *static_cast<superposition*>(pos->data);
			sp = { sp[rand() % sp.size()] };
		}
		void placeAndCollapse(gridState& gs, Grid::cell* pos, option opt) {
			superposition& sp = *static_cast<superposition*>(pos->data);
			sp = { (int)distance(gs.allOptions.begin(),find(gs.allOptions.begin(),gs.allOptions.end(),opt)) }; // todo
			propagate(gs, pos);
		}
	};

}

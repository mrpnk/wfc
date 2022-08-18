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
#include <filesystem>

#include "grid.h"

namespace wfc
{
	namespace fs = std::filesystem;

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


	// ---------------------------------------------------------------------------------------
	// The segments

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



	struct option {
		int modIdx;
		short rot;
		bool mirrorX;
		bool operator==(option const&)const = default;
	};


	class SegmentPalette{
		static const int m = 6; // tile size per axis in pixels
		static const int k = 5; // source texture size per axis in tiles
		static const bool MIRROR = 1;

		using seg = segment<m>;
		std::vector<seg> modus;
		std::vector<bool> fitlus;
		std::vector<option> allOptions;

		inline int faceElemIdx(int side, int i) const {
			int indices[4] = { i,m - 1 + i * m,m * m - 1 - i,(m - 1 - i) * m };
			return indices[side];
		}


		// Returns whether the two segments fit together with the given sides.
		bool fit(seg const& mod1, int s1, seg const& mod2, int s2, bool mirror) const {
			for (int i = 0; i < m; ++i)
				if (mod1.cols[faceElemIdx(s1, i)] != mod2.cols[faceElemIdx(s2, mirror ? i : m - 1 - i)])
					return false;
			return true;
		}

		// Computes the lookup table that stores which segments fit together.
		void computeFitLU() {
			AutoTimer at(g_timer, _FUNC_);
			for (int mod1 = 0; mod1 < modus.size(); ++mod1) {
				for (int s1 = 0; s1 < 4; ++s1) {
					for (int mod2 = 0; mod2 < modus.size(); ++mod2) {
						for (int s2 = 0; s2 < 4; ++s2) {
							for (int mirrorX = 0; mirrorX <= MIRROR; ++mirrorX)
								fitlus[(((mod1 * 4 + s1) * modus.size() + mod2) * 4 + s2) * (MIRROR + 1) + mirrorX] = fit(modus[mod1], s1, modus[mod2], s2, (bool)mirrorX);
						}
					}
				}
			}
		}

		inline bool fitlu(int mod1, int s1, int mod2, int s2, bool mirror) const {
			return fitlus[(((mod1 * 4 + s1) * modus.size() + mod2) * 4 + s2) * (MIRROR + 1) + mirror];
		}
		inline bool fitlu(int mod1, int rot1, bool mirror1, int side1, int mod2, int rot2, bool mirror2) const {
			int s1 = (side1 - rot1 + 4) % 4, s2 = (3 - side1 - rot2 + 4) % 4;
			if (mirror1) { s1 = (4 - s1) % 4; }
			if (mirror2) { s2 = (4 - s2) % 4; }
			return fitlu(mod1, s1, mod2, s2, mirror1 ^ mirror2);
		}

	public:
		void load(fs::path filename){
			AutoTimer at(g_timer, _FUNC_);
			const int nmodus = k * k;
			const int noptions = nmodus * 4 * (MIRROR + 1);
			modus.resize(nmodus);
			allOptions.resize(noptions);
			fitlus.resize((nmodus * 4) * (nmodus * 4) * (MIRROR + 1));

			if(fs::exists(filename)){
				std::ifstream infile(filename, std::ios::binary);
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
				for (int i = 0; i < modus.size(); ++i)
					for (short r = 0; r < 4; ++r)
						for (short mirr = 0; mirr <= MIRROR; ++mirr)
							allOptions[(i * 4 + r) * (MIRROR + 1) + mirr] = { i,r,(bool)mirr };
				computeFitLU();
			}
			else std::cout<< "Could not find segment file " << filename << std::endl;
		}
		option const& getOption(int i) const {
			return allOptions[i];
		}
		int getNumOptions() const {
			return allOptions.size();
		}
		inline bool isPossibleNeighbour(superposition const& sp, int dir, option const& other_opt) const {
			for (int i : sp)
				if (const option& own_opt = allOptions[i];
						fitlu(own_opt.modIdx, own_opt.rot, own_opt.mirrorX, dir,
						      other_opt.modIdx, other_opt.rot, other_opt.mirrorX))
					return true;
			return false;
		}
	};


	// ---------------------------------------------------------------------------------------
	// The possibilities

	struct gridState{

		Grid const* grid = nullptr;
		SegmentPalette const* palette = nullptr;
		wavefunction wave;
		superposition allSP;

		void initWithGrid(Grid const* g, SegmentPalette const* p){
			grid = g;
			palette = p;

			allSP.resize(palette->getNumOptions()); std::iota(std::begin(allSP), std::end(allSP), 0);

			// Create a superposition for each cell in the grid:
			wave.resize(grid->forAllCells([](auto) {}), allSP);
			grid->forAllCells([&, counter = 0](Grid::cell const& c)mutable{ c.data = &wave[counter++]; });
		}

		Grid::cell const* getPos(Grid::cell const* p, int d) const {
			return p->nn[d];
		}

		bool isCollapsed() const {
			for (auto& a : wave) if (a.size() == 0) return true; else if (a.size() > 1) return false;
			return true;
		}
		bool isStuck() const {
			for (auto& a : wave) if (a.size() == 0) return true;
			return false;
		}
		Grid::cell const* getLowestEntropyPos() const {
			std::vector<Grid::cell const*> minPosis; int minEntropy = std::numeric_limits<int>::max();
			grid->forAllCells([&](Grid::cell const& c) {
				superposition& sp = *static_cast<superposition*>(c.data);
				if (sp.size() > 1) {
					if (sp.size() < minEntropy) minPosis = { &c }, minEntropy = sp.size();
					else if (sp.size() == minEntropy) minPosis.push_back(&c);
				}
			});
			return minPosis[rand() % minPosis.size()];
		}
		bool isAnyMirrored() const {
			for (auto& a : wave)
				for (auto& b : a)
					if (palette->getOption(b).mirrorX)
						return true;
			return false;
		}


		void printCanvas(std::ostream& os, std::string sep = "") const {
			for (auto& sp : wave) {
				if (sp.size() != 1)
					os << "-1 0 0" << sep;
				else {
					const option& opt = palette->getOption(sp[0]);
					os << opt.modIdx << " " << opt.rot << " " << opt.mirrorX << sep;
				}
				os << std::endl;
			}
		}
	};



	// ---------------------------------------------------------------------------------------
	// The algorithm

	class WaveFunctionCollapser
	{
		gridState* state;
		SegmentPalette const* palette;

		void collapse(Grid::cell const* pos) {
			superposition& sp = *static_cast<superposition*>(pos->data);
			sp = { sp[rand() % sp.size()] };
		}
		void propagate(Grid::cell const* pos) {
			AutoTimer at(g_timer, _FUNC_);
			uniqueStack<Grid::cell const*> jobs;
			jobs.push(pos);
			while (!jobs.empty()) {
				Grid::cell const* p = jobs.pop();
				superposition& sp = *static_cast<superposition*>(p->data);
				for (int d = 0; d < 4; ++d) {
					if (Grid::cell const* nnp = state->getPos(p, d); nnp != nullptr) {
						superposition& nnsp = *static_cast<superposition*>(nnp->data);
						if (nnsp.size() <= 1) continue;
						if (std::erase_if(nnsp, [&](int& other_opt) {return !palette->isPossibleNeighbour(sp, d, palette->getOption(other_opt)); }))
							jobs.push(nnp);
					}
				}
			}
		}

	public:
		void solve(gridState* s) {
			AutoTimer at(g_timer, _FUNC_);

			state = s;
			palette = s->palette;

			while (!state->isCollapsed()) {
				auto p = state->getLowestEntropyPos();
				collapse(p);
				propagate(p);
			}

			std::cout << "isAnyMirrored = " << state->isAnyMirrored() << std::endl;
			std::cout << "isStuck = " << state->isStuck() << std::endl;
		}
	};

}

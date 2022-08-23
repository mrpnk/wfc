#pragma once

#include "grid.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include <vector>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <mutex>
#include <filesystem>
#include <random>


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
		int prio = 0;
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

		const Grid* grid = nullptr;
		SegmentPalette const* palette = nullptr;
		wavefunction wave;
		superposition allSP;

		void init(Grid const* g, SegmentPalette const* p){
			grid = g;
			palette = p;

			allSP.resize(palette->getNumOptions()); std::iota(std::begin(allSP), std::end(allSP), 0);

			reset();
		}

		// Creates the total superposition for each cell in the grid:
		void reset(){
			wave.resize(grid->forAllCells([](auto) {}), allSP);
			grid->forAllCells([&, counter = 0](face const& c)mutable{ c.data = &wave[counter++]; });
		}

		const face* getNeighbourFace(face const* p, int d) const {
			return p->nn[d]==-1?nullptr:&grid->faces[p->nn[d]];
		}

		bool isCollapsed() const {
			for (auto& a : wave) if (a.size() == 0) return true; else if (a.size() > 1) return false;
			return true;
		}
		bool isStuck() const {
			for (auto& a : wave) if (a.size() == 0) return true;
			return false;
		}
		const face& getLowestEntropyPos() const {
			std::vector<face const*> minPosis; int minEntropy = std::numeric_limits<int>::max();
			grid->forAllCells([&](face const& c) {
				superposition& sp = *static_cast<superposition*>(c.data);
				if (sp.size() > 1) {
					if (sp.size() < minEntropy) minPosis = { &c }, minEntropy = sp.size();
					else if (sp.size() == minEntropy) minPosis.push_back(&c);
				}
			});
			return *minPosis[rand() % minPosis.size()];
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
		std::mt19937 rg;

		// Fills the superposition of the slot based on hard conditions. The result is used as the initial state for wfc.
		void fillOptions(face const& c) {
			superposition& sp = *static_cast<superposition*>(c.data);
			sp = state->allSP;
		}

		void collapse(face const& c) {
			superposition& sp = *static_cast<superposition*>(c.data);
			sp = { sp[rand() % sp.size()] };
		}

		// Remove options of neighbouring cells that do not fit to any option of the given cell. Returns if the neighbours have options left.
		bool propagate(face const& c) {
			AutoTimer at(g_timer, _FUNC_);
			std::map<face const*, superposition> backup;
			uniqueStack<face const*> jobs;
			jobs.push(&c);
			while (!jobs.empty()) {
				face const* p = jobs.pop();
				superposition& sp = *static_cast<superposition*>(p->data);
				for (int d = 0; d < 4; ++d) {
					if (face const* nnp = state->getNeighbourFace(p, d); nnp != nullptr) {
						superposition& nnsp = *static_cast<superposition*>(nnp->data);
						if (nnsp.size() <= 1) continue;
						if (!backup.contains(nnp)) backup[nnp] = nnsp;

						erase_neighbours:
						if (std::erase_if(nnsp, [&](int& other_opt) {
							return !palette->isPossibleNeighbour(sp, d, palette->getOption(other_opt));
						}))
							jobs.push(nnp);

						if (nnsp.size() == 0) {
							if (!nnp->touched) {
								// we reached a dead end but there is hope: we can re-update this slot and we might get new options:
								//std::cerr << "re-update neighbour slot to save the wave..." << std::endl;
								fillOptions(*nnp);
								nnp->touched = true;
								if (nnsp != backup[nnp]) {
//									std::cerr << "re-update found new options: ";
//									for (auto& a: nnsp)
//										std::cerr << a << ",";
//									std::cerr << " ->Repeat." << std::endl;
									backup[nnp] = nnsp;
									goto erase_neighbours;
								}
							} else goto undo_propagate;
						}
					}
				}
			}
			return true;

			undo_propagate:
			for (auto& p: backup) {
				*static_cast<superposition*>(p.first->data) = p.second;
			}
			return false;
		}

		// Fills superpositions of all (touched) slots with the prefiltered options.
		void fillWave(bool all){
			state->grid->forAllCells([&, counter=0](face const& c)mutable{
				bool needsUpdate = c.touched||all;
				if(needsUpdate){// suppose this is a new slot
					fillOptions(c);
				}
			});
		}

		// Returns a list of numbers from 0 to n-1 in random order but such that the prio(ret[i]) >= prio(ret[i+1]).
		std::vector<int> sort_shuffle(int n, std::function<int(int)> prio, std::mt19937& rg){
			if(n==1)return {0};
			std::vector<int> rand_order(n);
			std::iota(std::begin(rand_order), std::end(rand_order), 0);
			std::shuffle(std::begin(rand_order), std::end(rand_order), rg);
			std::sort(std::begin(rand_order), std::end(rand_order), [&](int i, int j){return prio(i)>prio(j);});
			return rand_order;
		}

	public:
		void solve(gridState* s) {
			AutoTimer at(g_timer, _FUNC_);

			state = s;
			palette = s->palette;

			while (!state->isCollapsed()) {
				auto& p = state->getLowestEntropyPos();
				collapse(p);
				propagate(p);
			}

			std::cout << "isAnyMirrored = " << state->isAnyMirrored() << std::endl;
			std::cout << "isStuck = " << state->isStuck() << std::endl;
		}


		// Returns if successfull i.e. not stuck.
		bool backtrack_r(int& out_nrecursions, int& out_niters, int stage) {
			std::string ph = std::string(stage, '|');
			if (state->isCollapsed())
				return !state->isStuck();
			out_nrecursions++;
			const auto& c = state->getLowestEntropyPos();
			superposition& sp = *static_cast<superposition*>(c.data);
			//	std::cerr<<"getLowestEntropyPos = "<<p->mycell<< " ori "<<p->ori<<endl;
			if (sp.size() == 1)
				c.touched = false;
			superposition backup = sp;
			//std::cerr << ph << "backup="; for (int i : backup) std::cerr << i << ", ";
			//std::cerr << std::endl;
			std::vector<int> rand_order = sort_shuffle(backup.size(),
													   [&](int i) {return palette->getOption(i).prio; }, rg);
			for (int idx : rand_order) {
				// Make a choice:
				int a = backup[idx];
				sp = { a };
				//std::cerr << ph << "choose " << a << std::endl;
				if (!propagate(c))
				{
					//std::cerr << ph << "redone choice " << a << std::endl;
					continue;
				}
				//std::cerr << ph << "recursion ->" << std::endl;
				if (!backtrack(out_nrecursions, out_niters, stage + 1))
				{
					//std::cerr << ph << "backtrace " << a << std::endl;
					continue;
				}
				//std::cerr << ph << "success" << std::endl;
				return true;
			}
			sp = backup;
			return false;
		}

		bool backtrack(int& out_nrecursions, int& out_niters, int stage) {

			std::stack<int> jobs;
			jobs.push(0);
			while(jobs.empty()){
				jobs.pop();

				std::string ph = std::string(stage, '|');
				if (state->isCollapsed())
					return !state->isStuck();
				out_nrecursions++;
				const auto& c = state->getLowestEntropyPos();
				superposition& sp = *static_cast<superposition*>(c.data);
				//	std::cerr<<"getLowestEntropyPos = "<<p->mycell<< " ori "<<p->ori<<endl;
				if (sp.size() == 1)
					c.touched = false;
				superposition backup = sp;
				//std::cerr << ph << "backup="; for (int i : backup) std::cerr << i << ", ";
				//std::cerr << std::endl;
				std::vector<int> rand_order = sort_shuffle(backup.size(),
				                                           [&](int i) {return palette->getOption(i).prio; }, rg);
				for (int idx : rand_order) {
					// Make a choice:
					int a = backup[idx];
					sp = { a };
					//std::cerr << ph << "choose " << a << std::endl;
					if (!propagate(c))
					{
						//std::cerr << ph << "redone choice " << a << std::endl;
						continue;
					}
					//std::cerr << ph << "recursion ->" << std::endl;

					if (!backtrack(out_nrecursions, out_niters, stage + 1))
					{
						//std::cerr << ph << "backtrace " << a << std::endl;
						continue;
					}
					//std::cerr << ph << "success" << std::endl;
					return true;
				}
				sp = backup;
				return false;
			}
		}


		// Solves the wave with backtracking.
		void solveRecursive(gridState* s) { //, newModuleCB_t&& newMod
			AutoTimer at(g_timer, _FUNC_);

			state = s;
			palette = s->palette;

			fillWave(true);
			auto action = [&]() mutable {

//				// Map the intern wave function to the grid slots:
//				wavefunction gridwave(gs..getNumCells() * 8);
//				g.forAllCells([&, counter = 0](Cell&, int, slot& s)mutable{
//					gridwave[counter++] = &s.sp;
//				});

				// Collapse and propagate:
				auto t0 = std::chrono::high_resolution_clock::now();
				int nrecursions = 0, niters = 0;

				backtrack(nrecursions, niters, 0);

				std::cout << "Wave function collapse done:" << std::endl;
				std::cout << "> Time = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << " ms" << std::endl;
				std::cout << "> Wave size = " << state->wave.size() << std::endl;
				std::cout << "> recursion = " << nrecursions << std::endl;
				std::cout << "> neighbour-interactions = " << niters << std::endl;
				std::cout << "> isAnyMirrored = " << state->isAnyMirrored() << std::endl;
				std::cout << "> isStuck = " << state->isStuck() << std::endl;

				//realize(g, newMod);
			};

			action();
		}
	};

}

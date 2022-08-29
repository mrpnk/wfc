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
#include <execution>

namespace wfc
{
	template<typename T>
	class uniqueStack {
		std::stack<T> st;
		std::set<T> se;
		mutable std::mutex mtx;
	public:
		bool push(T i) {
			//AutoTimer at(g_timer, _FUNC_);
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

	using seg_id = int;
	using trafo_id = char;
	using edge_id = int;


	struct transform2D{
		const int N;
		const bool reflectX; // first reflect, then rotate
		const int rot;
		const edge_id* const mym;

		static inline edge_id m3[3*2*6], m4[4*2*8], m5[5*2*10], m6[6*2*12];
		static constexpr edge_id (*ms[7])[] = {0,0,0,&m3,&m4,&m5,&m6};

		transform2D(int n, bool refl, int rot)
			: N{n}, reflectX{refl}, rot{rot}, mym{(*ms[N])+(rot*2+reflectX)*(2*N)}{
			assert(3<=n && n<=6);
		}
		inline edge_id apply(edge_id e) const {
			// With careless implementation, a single integer multiplication in this line
			// can be responsible for ~10% of total execution time! The runtime-N really hurts here.
			// I think, this solution is near optimal, though.
			return mym[e];
		}
		static void generate(){
			for(int N = 3; N <= 6; ++N){
				for(edge_id e = 0; e < 2*N; ++e){
					for(int m = 0; m < 2; ++m){
						for(int r = 0; r < N; ++r){
							int off = e/N;
							edge_id x = (e%N-r+N)%N+off*N;
							if(m) x = ((N - x%N) % N)+(1-off)*N;
							(*ms[N])[(r*2+m)*(2*N)+e] = x;
						}
					}
				}
			}
		}
	};

	/*          cell/slot
	 *             /\
	 *            /  \
	 *           / /\ \
	 *      e=1 / /  \ \ e=2
	 *         / /    \ \
	 *        / / seg  \ \
	 *       / /________\ \
	 *      /______________\
	 *            e=0
	 */
	struct symmetryGroup2D {
		virtual wfc::trafo_id getNumElems() const = 0;
		virtual wfc::transform2D getTrafo(wfc::trafo_id) const = 0;
	};
	template<int N, bool reflect = true>
	struct dihedralGroup : public symmetryGroup2D {
		wfc::trafo_id getNumElems() const override { return N*(1+reflect); }
		wfc::transform2D getTrafo(wfc::trafo_id i) const override{
			return wfc::transform2D(N,i / N > 0,i % N);
		}
	};


	struct option {
		seg_id segment;
		transform2D trafo;
		int prio = 0;
		bool operator==(option const&)const = default;
	};

	using option_id = short;
	using superposition = std::vector<option_id>;
	using wavefunction = std::vector<superposition>;

	class SegmentPalette {
	protected: // These members are fixed to avoid vftable lookups for better performance.
		std::vector<wfc::option> options;
		std::vector<bool> fitlus;
		int nSegs = 0, nTrafos = 0; // Need to be set by the derived class

		inline size_t luidx(wfc::seg_id s1, wfc::edge_id e1, wfc::seg_id s2, wfc::edge_id e2) const {
			return ((s1 * nTrafos + e1) * nSegs + s2) * nTrafos + e2;
		}

		virtual bool fit(seg_id s1, edge_id e1, seg_id s2, edge_id e2) const = 0;

		// Computes the lookup table that stores which segments fit together.
		void computeFitLU() {
			AutoTimer at(g_timer, _FUNC_);
			fitlus.resize(std::pow(nSegs*nTrafos, 2)); // This is too many TODO /2
			for (option const& o1 : options) {
				for (option const& o2 : options) {
					wfc::edge_id e1 = o1.trafo.apply(0);
					wfc::edge_id e2 = o2.trafo.apply(0);
					fitlus[luidx(o1.segment, e1, o2.segment, e2)] = fit(o1.segment, e1, o2.segment, e2);
				}
			}
		}


	public:
		int getNumOptions() const {
			return options.size();
		}
		inline option const& getOption(wfc::option_id i) const {
			return options[i];
		}

		// Fills in the options for n-gons.
		void getOptions(std::vector<option_id>& opts, int nEdges) const {
			for(option_id id = 0; option const& o : options){
				if(o.trafo.N == nEdges)
					opts.push_back(id);
				++id;
			}
		}


		// Returns whether o2 fits to any option in sp. This is performance-critical!
		inline bool isPossibleNeighbour(wfc::superposition const& sp, edge_id e1, wfc::option const& o2, edge_id e2) const {
			e2 = o2.trafo.apply(e2);
			for (option_id i : sp){
				const wfc::option& o1 = options[i];
				if (fitlus[luidx(o1.segment,o1.trafo.apply(e1), o2.segment,e2)])
					return true;
			}
			return false;
		}
	};

	// ---------------------------------------------------------------------------------------
	// The possibilities

	template<class grid_t>
	class gridState{
	public:
		using face_t = typename grid_t::face_t;
		const grid_t* grid = nullptr;
		const SegmentPalette* palette = nullptr;

		wavefunction wave;
		std::vector<superposition> allSPn; // one for each edge count

		struct faceInfo{
			bool dirty = false;
			bool disabled = false;
		};
		std::vector<faceInfo> faceInfos;


		void init(grid_t const* g, SegmentPalette const* p){
			grid = g;
			palette = p;

			// Creates the total superpositions for each edge count
			allSPn.resize(7);
			for(int n = 3; n <= 6; ++n)
				palette->getOptions(allSPn[n], n);

			reset();
		}

		void reset(){
			// Allow each face all options for its edge count
			wave.resize(grid->faces.size());
			faceInfos.resize(grid->faces.size());
			for(int i = 0; i < grid->faces.size(); ++i) {
				int nEdges = grid->faces[i].getNumCorners();
				wave[i] = allSPn[nEdges];
				faceInfos[i].dirty = false;
				faceInfos[i].disabled = (nEdges < 3);
			}
		}

		// Fills the superposition of the slot based on hard conditions, e.g. by marching cubes.
		// The result is used as the initial state for wfc.
		void fillOptions(face_id f) {
			getOptions(f) = allSPn[grid->faces[f].getNumCorners()];
		}

		// Fills superpositions of all (touched) slots with the prefiltered options.
		void fillWave(bool all){
			grid->forAllFaces([&, counter=0](face_t const& c)mutable{
				bool needsUpdate = getInfo(&c).dirty || all;
				if(needsUpdate){// suppose this is a new slot
					fillOptions(c);
				}
			});
		}

		face_id getLowestEntropyFace() const {
			AutoTimer at(g_timer, _FUNC_);

//			int i;
//			do{(i = rand() % grid->faces.size());}
//			while(this->getOptions(i).size() <= 1);
//			return i;

			std::vector<face_id> minPosis; int minEntropy = std::numeric_limits<int>::max();
			grid->forAllFaces([&, this](face_t const& c) {
				face_id id = &c-grid->faces.data();
				const superposition& sp = this->getOptions(id);
				if (sp.size() > 1) {
					if (sp.size() < minEntropy) minPosis = { id }, minEntropy = sp.size();
					else if (sp.size() == minEntropy) minPosis.push_back(id);
				}
			});
			return minPosis[rand() % minPosis.size()];
		}

		const face_id getNeighbourFace(face_id f, int d) const {
			return grid->faces[f].neighbours[d];
		}
		auto& getOptions(this auto&& self, face_id f) {
			return self.wave[f];
		}

		auto& getInfo(this auto&& self, face_id f) {
			return self.faceInfos[f];
		}

		bool isCollapsed() const {
			for (int i = 0; auto& a : wave) {
				if(!faceInfos[i].disabled){if (a.size() == 0) return true; else if (a.size() > 1) return false;}
				++i;
			}
			return true;
		}
		bool isStuck() const {
			for (int i = 0; auto& a : wave){
				if(!faceInfos[i].disabled){if (a.size() == 0) return true;}
				++i;
			}
			return false;
		}
		bool isAnyReflected() const {
			for (auto& a : wave)
				for (auto& b : a)
					if (palette->getOption(b).trafo.reflectX)
						return true;
			return false;
		}

		void checkSolution() const{
			assert(isCollapsed() && !isStuck());
			int nDefects = 0;
			for(int f0 = 0; f0 < grid->faces.size(); ++f0){
				if(grid->faces[f0].getNumCorners() < 3) continue;
				for(int d = 0; d < grid->faces[f0].getNumCorners(); ++d){
					int f1 = grid->faces[f0].neighbours[d];
					if(f1 < 0) continue;
					if(grid->faces[f1].getNumCorners() < 3) continue;
					bool fits = palette->isPossibleNeighbour(wave[f0], d, palette->getOption(wave[f1].front()), grid->getComplFace(f0, d));
					if(!fits) nDefects++;
				}
			}
			std::cerr << "Num fit Defects = "<< nDefects <<std::endl;
		}

		void print(std::ostream& os, std::string sep = "") const {
			for (auto& sp : wave) {
				if (sp.size() != 1)
					os << "-1 0 0" << sep;
				else {
					const option& opt = palette->getOption(sp[0]);
					os << opt.segment << " " << opt.trafo.rot << " " << (int)opt.trafo.reflectX << sep;
				}
				os << std::endl;
			}
		}
	};



	// ---------------------------------------------------------------------------------------
	// The algorithm

	template<class grid_t>
	class WaveFunctionCollapser
	{
		using face_t = typename grid_t::face_t;
		using state_t = typename gridState<grid_t>;
		state_t* state;
		SegmentPalette const* palette;
		std::mt19937 rg;

		// Returns a list of numbers from 0 to n-1 in random order but such that the prio(ret[i]) >= prio(ret[i+1]).
		std::vector<int> sort_shuffle(int n, std::function<int(int)> prio, std::mt19937& rg){
			if(n==1)return {0};
			std::vector<int> rand_order(n);
			std::iota(std::begin(rand_order), std::end(rand_order), 0);
			std::shuffle(std::begin(rand_order), std::end(rand_order), rg);
			std::sort(std::begin(rand_order), std::end(rand_order), [&](int i, int j){return prio(i)>prio(j);});
			return rand_order;
		}


		void collapse(face_id f) {
			superposition& sp = state->getOptions(f);
			sp = { sp[rand() % sp.size()] };
			//std::cerr << "Collapsed "<<f<<" to "<<sp.front()<<std::endl;
		}

		// Remove options of neighbouring cells that do not fit to any option of the given cell. Returns if the neighbours have options left.
		bool propagate(face_id f) {
			AutoTimer at(g_timer, _FUNC_);

			//auto hardBackup = state->wave;

			std::map<face_id, superposition> backup;
			uniqueStack<face_id> jobs;
			const auto& faces = state->grid->faces;
			jobs.push(f);
			while (!jobs.empty()) {
				face_id f0 = jobs.pop();
				const superposition& sp0 = state->getOptions(f0);
				for (edge_id d0 = 0; d0 < faces[f0].getNumCorners(); ++d0) {
					face_id f1 = state->getNeighbourFace(f0, d0);
					if (f1 < 0) continue;

					superposition& sp1 = state->getOptions(f1);
					if (sp1.size() <= 1) continue;
					backup.try_emplace(f1, sp1);

					erase_neighbours:
					{
						//AutoTimer at(g_timer, "erase neighbours");
						const int d1 = state->grid->getComplFace(f0, d0);
						if (std::erase_if(sp1, [&](option_id const& o1) {
							return !palette->isPossibleNeighbour(sp0, d0, palette->getOption(o1), d1);
						}))
							jobs.push(f1);
					}

					if (sp1.empty()) {
						std::cerr << "empty" << std::endl;

						if (auto& fi = state->getInfo(f1); !fi.dirty) {
							// we reached a dead end but there is hope: we can re-update this slot and we might get new options:
							state->fillOptions(f1);
							//fi.dirty = true;
							if (sp1 != backup[f1]) {
								std::cerr << "saved: found new options" << std::endl;
								backup.insert_or_assign(f1, sp1);
								goto erase_neighbours;
							}
							else
								goto undo_propagate;
						} else goto undo_propagate;
					}
				}
			}
			return true;

			undo_propagate:
			std::cerr<< "undo propagate: " <<backup.size() << std::endl;
			{
				AutoTimer at(g_timer, "undo propagate");
				for (const auto& p : backup)
					state->getOptions(p.first) = p.second;
				//assert(state->wave == hardBackup);
			}
			return false;
		}

		// Returns if successful i.e. not stuck.
		bool backtrack(int& out_nrecursions, int& out_niters, int stage) {
			std::string ph = std::string(stage, '|');
			if (state->isCollapsed())
				return !state->isStuck();
			out_nrecursions++;
			const auto& c = state->getLowestEntropyFace();
			superposition& sp = state->getOptions(&c);
			//	std::cerr<<"getLowestEntropyPos = "<<p->mycell<< " ori "<<p->ori<<endl;
			if (sp.size() == 1)
				state->getInfo(&c).dirty = false;
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

	public:
		void solve(state_t* s) {
			AutoTimer at(g_timer, _FUNC_);

			state = s;
			palette = s->palette;
			superposition backup;

			while (!state->isCollapsed()) {
				face_id f = state->getLowestEntropyFace();
				backup = state->getOptions(f);
				collapse(f);
				if(!propagate(f))
					state->getOptions(f) = backup;
			}

			std::cout << "isAnyReflected = " << state->isAnyReflected() << std::endl;
			std::cout << "isStuck = " << state->isStuck() << std::endl;
		}

		// Solves the wave with backtracking.
		void solveRecursive(state_t* s) {
			AutoTimer at(g_timer, _FUNC_);

			state = s;
			palette = s->palette;

			state->fillWave(true);
			auto action = [&]() mutable {

//				// Map the intern wave function to the grid slots:
//				wavefunction gridwave(gs..getNumCells() * 8);
//				g.forAllFaces([&, counter = 0](Cell&, int, slot& s)mutable{
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
				std::cout << "> isAnyReflected = " << state->isAnyReflected() << std::endl;
				std::cout << "> isStuck = " << state->isStuck() << std::endl;

				//realize(g, newMod);
			};

			action();
		}
	};

}

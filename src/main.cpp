#include "timer.h"
#include "grid.h"
#include "wfc.h"

#include <regex>

struct colour {
	float r, g, b;
	bool operator<(colour const& c) const { return r < c.r ? true : r > c.r ? false : g < c.g ? true : g > c.g ? false : b < c.b; }
	bool operator==(colour const& c) const = default;
	friend std::ostream& operator<<(std::ostream& os, colour const& col) { return os << col.r << " " << col.g << " " << col.b; }
};


class ImagePalette : public wfc::SegmentPalette{
	// first reflect, then rotate
	static const int m = 6; // tile size per axis in pixels
	static const int k = 5; // source texture size per axis in tiles
	static const bool MIRROR = 1;

	template<int m> struct segment {
		colour cols[m * m];
	};

	using seg = segment<m>;
	std::vector<seg> modus;
	std::vector<bool> fitlus;
	std::vector<wfc::option> allOptions;

	/*           side=0
	 *           ----->
	 *         A 0 1 2 |
	 *  side=3 | 3 4 5 | side=1
	 *         | 6 7 8 V
	 *           <-----
	 *           side=2
	 */
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
						for (int mirrorX = 0; mirrorX <= (short)MIRROR; ++mirrorX)
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
	void load(std::filesystem::path filename){
		AutoTimer at(g_timer, _FUNC_);
		const int nmodus = k * k;
		const int noptions = nmodus * 4 * (MIRROR + 1);
		modus.resize(nmodus);
		allOptions.resize(noptions);
		fitlus.resize((nmodus * 4) * (nmodus * 4) * (MIRROR + 1));

		if(std::filesystem::exists(filename)){
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
					for (short mirr = 0; mirr <= (short)MIRROR; ++mirr)
						allOptions[(i * 4 + r) * (MIRROR + 1) + mirr] = { i,r,(bool)mirr };
			computeFitLU();
		}
		else std::cout<< "Could not find segment file " << filename << std::endl;
	}
	int getNumOptions() const override {
		return allOptions.size();
	}
	wfc::option const& getOption(wfc::option_id i) const override {
		return allOptions[i];
	}
	bool isPossibleNeighbour(wfc::superposition const& sp, int dir, wfc::option const& otherOpt, int otherDir) const override {
		for (int i : sp)
			if (const wfc::option& own_opt = allOptions[i];
					fitlu(own_opt.modIdx, own_opt.rot, own_opt.mirrorX, dir,
					      otherOpt.modIdx, otherOpt.rot, otherOpt.mirrorX))
				return true;
		return false;
	}
};

namespace fs = std::filesystem;

class CarcassonnePalette : public wfc::SegmentPalette{

	using interface_id = unsigned char;

	struct seg {
		interface_id interfaces[4];
		bool reflectable;
	};

	std::vector<seg> segments;
	std::vector<bool> fitlus;
	std::vector<wfc::option> allOptions;
	const int MIRROR = 0;

	// Returns whether the two segments fit together with the given sides.
	bool fit(seg const& s1, int e1, seg const& s2, int e2, bool mirror) const {
//		e1 = (e1+1)%4;
//		e2 = (e2+1)%4;
		return s1.interfaces[e1] == s2.interfaces[e2];
	}

	// Computes the lookup table that stores which segments fit together.
	void computeFitLU() {
		AutoTimer at(g_timer, _FUNC_);
		for (int s1 = 0; s1 < segments.size(); ++s1) {
			for (int e1 = 0; e1 < 4; ++e1) {
				for (int s2 = 0; s2 < segments.size(); ++s2) {
					for (int e2 = 0; e2 < 4; ++e2) {
						for (int mirrorX = 0; mirrorX <= (short)MIRROR; ++mirrorX)
							fitlus[(((s1 * 4 + e1) * segments.size() + s2) * 4 + e2) * (MIRROR + 1) + mirrorX] = fit(segments[s1], e1, segments[s2], e2, (bool)mirrorX);
					}
				}
			}
		}
	}

	inline bool fitlu(int mod1, int s1, int mod2, int s2, bool mirror) const {
		return fitlus[(((mod1 * 4 + s1) * segments.size() + mod2) * 4 + s2) * (MIRROR + 1) + mirror];
	}
	inline bool fitlu(int seg1, int rot1, bool refl1, int side1,
	                  int seg2, int rot2, bool refl2, int side2) const {
		int e1 = (side1 - rot1 + 4) % 4, e2 = (side2 - rot2 + 4) % 4;
		if (refl1) { e1 = (4 - e1) % 4; }
		if (refl2) { e2 = (4 - e2) % 4; }
		return fitlu(seg1, e1, seg2, e2, refl1 ^ refl2);
	}

public:
	void load(fs::path directory){
		AutoTimer at(g_timer, _FUNC_);

		const auto mask = std::regex("[012]{4}(_\\d)?\\.png");
		int nSegments = 0, nOptions = 0;
		for (const auto& file : fs::directory_iterator(directory)){
			auto filename = file.path().filename().string();
			if(file.is_regular_file() && std::regex_match(filename, mask)){
				auto& seg = segments.emplace_back();
				for(int i = 0; i < 4; ++i) seg.interfaces[i] = filename[i];
				seg.reflectable = false;
				nSegments += 1;
				nOptions += 4; // todo!

				std::cout << "Found segment file " << filename << std::endl;
			}
		}

		allOptions.resize(nOptions);
		fitlus.resize(nOptions * nOptions);

		for (int i = 0; i < segments.size(); ++i)
			for (short r = 0; r < 4; ++r)
				for (short mirr = 0; mirr <= (short)MIRROR; ++mirr)
					allOptions[(i * 4 + r) * (MIRROR + 1) + mirr] = { i,r,(bool)mirr };
		computeFitLU();
	}
	int getNumOptions() const override {
		return allOptions.size();
	}
	const wfc::option& getOption(wfc::option_id i) const override {
		return allOptions[i];
	}
	inline bool isPossibleNeighbour(wfc::superposition const& sp, int dir, wfc::option const& otherOpt, int otherDir) const override {
		for (wfc::option_id i : sp)
			if (const wfc::option& own_opt = allOptions[i];
					fitlu(own_opt.modIdx, own_opt.rot, own_opt.mirrorX, dir,
					      otherOpt.modIdx, otherOpt.rot, otherOpt.mirrorX, otherDir))
				return true;
		return false;
	}
};


int main()
{
	srand(273);

	// In the beginning there was the grid
	using grid_t = grid<upto<6>, always<4>>;
	grid_t primal;

	HexQuadGenerator generator;
	generator.generate(50);
	generator.relax(10);

	generator.print();

	generator.convert(primal);

//	using grid_t = grid<always<4>, always<4>>;
//	grid_t primal;
//
//	SquareGenerator generator;
//	generator.generate(30);
//	generator.convert(primal);

	std::ofstream file("grid.txt");
	primal.print(file);
	file.close();

	// Compute the dual because it looks so nice
	grid_t::dual_t dual;
	computeDualGrid(primal, dual);

	std::cout << "primal defects: " << std::endl;
	checkIntegrity(primal);

	std::cout << "dual defects: " << std::endl;
	checkIntegrity(dual);

	std::ofstream dualfile("dual.txt");
	dual.print(dualfile);
	dualfile.close();


	// Load the segments from file
	ImagePalette palette;
	palette.load("modules.txt");
//	CarcassonnePalette palette;
//	palette.load("carcassonne/", primal.meta);

	// Create the option space for this grid
	wfc::gridState<grid_t> state;
	state.init(&primal, &palette);


	// Find a constellation that satisfies all conditions
	wfc::WaveFunctionCollapser<grid_t> collapser;
	collapser.solve(&state);
	//collapser.solveRecursive(&state);

	std::ofstream outfile("canvas.txt");
	state.printCanvas(outfile, " ");
	outfile.close();

	g_timer.print();
}

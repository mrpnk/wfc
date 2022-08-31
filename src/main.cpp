#include "timer.h"
#include "grid.h"
#include "wfc.h"

#include <regex>

//#define BENCHMARK


class ImagePalette : public wfc::SegmentPalette{
	static const int m = 6; // tile size per axis in pixels
	static const int k = 5; // source texture size per axis in tiles
	static const bool allowReflections = true;
	struct colour {
		float r, g, b;
		bool operator<(colour const& c) const { return r < c.r ? true : r > c.r ? false : g < c.g ? true : g > c.g ? false : b < c.b; }
		bool operator==(colour const& c) const = default;
		friend std::ostream& operator<<(std::ostream& os, colour const& col) { return os << col.r << " " << col.g << " " << col.b; }
	};

	template<int m> struct segment {
		colour cols[m * m];
		const wfc::symmetryGroup2D* symm;

		/*           side=0
		 *           ----->
		 *         A 0 1 2 |
		 *  side=3 | 3 4 5 | side=1
		 *         | 6 7 8 V
		 *           <-----
		 *           side=2
		 */
		inline colour faceElem(wfc::edge_id e, int i) const {
			int indices[8] = { i, m-1 + i*m, m*m-1-i, (m-1-i)*m, // CW
			                   m-1-i, m*m-1-i*m, m*m-m+i, i*m};  // CCW
			return cols[indices[e]];
		}
	};
	using seg = segment<m>;
	std::vector<seg> segments;

	const wfc::dihedralGroup<4, allowReflections> d4;


	// Returns whether the two segments fit together with the given edges.
	bool fit(wfc::seg_id s1, wfc::edge_id e1, wfc::seg_id s2, wfc::edge_id e2) const override {
		for (int i = 0; i < m; ++i)
			if (segments[s1].faceElem(e1, i) != segments[s2].faceElem(e2, m-1-i))
				return false;
		return true;
	}

public:
	void load(std::filesystem::path filename){
		AutoTimer at(g_timer, _FUNC_);
		nSegs = k * k;
		nTrafos = d4.getNumElems();
		segments.resize(nSegs);

		if(std::filesystem::exists(filename)){
			colour imagedata[m * m * k * k];
			std::ifstream infile(filename, std::ios::binary);
			infile.read((char*)(imagedata), sizeof(imagedata));
			infile.close();
			for (int i = 0; i < k; ++i) {
				for (int j = 0; j < k; ++j) {
					for (int a = 0; a < m; ++a)
						for (int b = 0; b < m; ++b)
							segments[(i * k + j)].cols[a * m + b] = imagedata[m * m * k * i + m * j + m * k * a + b];
				}
			}
			for (wfc::seg_id s = 0; s < segments.size(); ++s) {
				auto& sy = segments[s].symm = &d4;
				for (wfc::trafo_id tr = 0; tr < sy->getNumElems(); ++tr)
					options.push_back({.segment=s, .trafo=sy->getTrafo(tr), .prio=0});
			}
			computeFitLU();
		}
		else std::cout << "Could not find segment file " << filename << std::endl;
	}
};


class CarcassonnePalette : public wfc::SegmentPalette{
	using interface_id = unsigned char;
	static const bool allowReflections = false;

	struct seg {
		interface_id interfaces[12]; // max 6 edges and two directions each
		const wfc::symmetryGroup2D* symm;
	};
	std::vector<seg> segments;

	const wfc::dihedralGroup<3, allowReflections> d3;
	const wfc::dihedralGroup<4, allowReflections> d4;
	const wfc::dihedralGroup<5, allowReflections> d5;
	const wfc::dihedralGroup<6, allowReflections> d6;

	// Returns whether the two segments fit together with the given sides.
	bool fit(wfc::seg_id s1, wfc::edge_id e1, wfc::seg_id s2, wfc::edge_id e2) const override {
		assert(e1 < segments[s1].symm->getNumElems());
		assert(e2 < segments[s2].symm->getNumElems());
		return segments[s1].interfaces[e1] == segments[s2].interfaces[e2];
	}

public:
	void load(std::filesystem::path directory){
		AutoTimer at(g_timer, _FUNC_);
		namespace fs = std::filesystem;

		nSegs = 0;
		nTrafos = d6.getNumElems(); // use the maximum

		std::vector<std::string> filenames;
		const auto mask = std::regex("[012]{3,6}(_\\d)?\\.png");
		for (const auto& file : fs::directory_iterator(directory)) {
			auto filename = file.path().filename().string();
			if (file.is_regular_file() && std::regex_match(filename, mask)) {
				filenames.push_back(filename);
			}
		}

		std::ofstream file(directory/"carcaIndex.txt");
		for(auto& a : filenames) file << a << "\n";
		file.close();

		for(auto filename : filenames){
			std::cout << "Found segment file " << filename << std::endl;
			auto& seg = segments.emplace_back();
			wfc::seg_id id = segments.size()-1;
			int nEdges = filename.find_first_of("_.");
			for(wfc::edge_id i = 0; i < nEdges; ++i) seg.interfaces[i] = seg.interfaces[nEdges+i] = filename[i];

			// Assign a symmetry group and generate options
			switch (nEdges) {
				case 3: seg.symm = &d3; break;
				case 4: seg.symm = &d4; break;
				case 5: seg.symm = &d5; break;
				case 6: seg.symm = &d6; break;
			}
			for (wfc::trafo_id tr = 0; tr < seg.symm->getNumElems(); ++tr)
				options.push_back({.segment=id, .trafo=seg.symm->getTrafo(tr), .prio=0});
			nSegs += 1;
		}

		computeFitLU();
	}
};



int main()
{
	srand(273);
	wfc::transform2D::generate();

#if defined BENCHMARK
	{
		// In the beginning there was the grid
		using grid_t = grid<upto<6>, always<4>>;
		grid_t primal;

		HexQuadGenerator generator;
		generator.generate(50);
		generator.convert(primal);

		std::ofstream file("grid.txt");
		primal.print(file);
		file.close();


		// Load the segments from file
		ImagePalette palette;
		palette.load("modules.txt");


		// Create the option space for this grid
		wfc::gridState<grid_t> state;
		state.init(&primal, &palette);

		// Find a constellation that satisfies all conditions
		wfc::WaveFunctionCollapser<grid_t> collapser;
		collapser.solve(&state);


		std::ofstream outfile("canvas.txt");
		state.print(outfile);
		outfile.close();

		g_timer.print();

		return 0;
	};
#endif

	using grid_t = grid<upto<6>, always<4>>;
	grid_t primal;

	HexQuadGenerator generator;
	generator.generate(4);
	generator.relax(10);
	generator.convert(primal);


	std::ofstream file("grid.txt");
	primal.print(file);
	file.close();

	std::ofstream debugfile("gridDebug.txt");
	primal.printDebug(debugfile);
	debugfile.close();



	// Compute the dual because it looks so nice
	grid_t::dual_t dual;
	computeDualGrid(primal, dual);

	std::cout << "primal defects: ";
	checkIntegrity(primal);

	std::cout << "dual defects: ";
	checkIntegrity(dual);

	std::ofstream dualfile("dual.txt");
	dual.print(dualfile);
	dualfile.close();

	std::ofstream dualDebugfile("dualDebug.txt");
	dual.printDebug(dualDebugfile);
	dualDebugfile.close();


	// Load the segments from file
//	ImagePalette palette;
//	palette.load("modules.txt");
	CarcassonnePalette palette;
	palette.load("carcassonne/");



//	// Create the option space for this grid
//	wfc::gridState<grid_t> state;
//	state.init(&primal, &palette);
//
//	// Find a constellation that satisfies all conditions
//	wfc::WaveFunctionCollapser<grid_t> collapser;
//	collapser.solve(&state);



	wfc::gridState<grid_t::dual_t> state;
	state.init(&dual, &palette);

	// Find a constellation that satisfies all conditions
	wfc::WaveFunctionCollapser<grid_t::dual_t> collapser;
	collapser.solve(&state);
	state.checkSolution();


	std::ofstream outfile("canvas.txt");
	state.print(outfile);
	outfile.close();

	g_timer.print();
}

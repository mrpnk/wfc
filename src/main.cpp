#include "grid.h"
#include "wfc.h"

int main()
{
	srand(27);

	// In the beginning there was the grid
	Grid g;
	g.init(10, 5);

	std::ofstream file("grid.txt");
	g.print(file);
	file.close();

	// Load the segments from file
	wfc::gridState gs;
	gs.load();
	gs.initWithGrid(g);

	wfc::WaveFunctionCollapser collapser;
	collapser.solve(g, gs);

	std::ofstream outfile("canvas.txt");
	gs.printCanvas(outfile, " ");
	outfile.close();
}

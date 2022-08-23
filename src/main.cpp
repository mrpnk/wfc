#include "timer.h"
#include "grid.h"
#include "wfc.h"

int main()
{
	srand(273);

	// In the beginning there was the grid
	Grid grid;
	GridGenerator gg;
	gg.construct(50);
	gg.relax(10);
	gg.convert(grid);

	std::ofstream file("grid.txt");
	grid.print(file);
	file.close();


	// Load the segments from file
	wfc::SegmentPalette palette;
	palette.load("modules.txt");


	// Create the option space for this grid
	wfc::gridState state;
	state.init(&grid, &palette);


	// Find a constellation that satisfies all conditions
	wfc::WaveFunctionCollapser collapser;
	collapser.solve(&state);
	//collapser.solveRecursive(&state);

	std::ofstream outfile("canvas.txt");
	state.printCanvas(outfile, " ");
	outfile.close();

	g_timer.print();
}

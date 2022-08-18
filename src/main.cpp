#include "timer.h"
#include "grid.h"
#include "wfc.h"

int main()
{
	srand(23);

	// In the beginning there was the grid
	Grid grid;
	grid.init(5);
	grid.relax(5);

	std::ofstream file("grid.txt");
	grid.print(file);
	file.close();


	// Load the segments from file
	wfc::SegmentPalette palette;
	palette.load("modules.txt");


	// Create the option space for this grid
	wfc::gridState state;
	state.initWithGrid(&grid, &palette);


	// Find a constellation that satisfies all conditions
	wfc::WaveFunctionCollapser collapser;
	collapser.solve(&state);

	std::ofstream outfile("canvas.txt");
	state.printCanvas(outfile, " ");
	outfile.close();

	g_timer.print();
}

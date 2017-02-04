#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <Utils/FDMGridMPI.h>

#include "TaskFunctions.h"
#include "Solvers/Functions.h"

int main()
{
	size_t nGridParam, nBWidth;

	// set the partition number of the square side
	nGridParam = 4;
	// set maximal number of nonzero matrix elements in row
	nBWidth = 5;

	// create grid
	fdm_grid grid(nGridParam);

	// allocate SLAE
	// - sparse matrix
	sparse_matrix mA(grid.nodes_number(), nBWidth);
	// - solution and right part
	vector vSolution(grid.nodes_number()), vRightPart(grid.nodes_number());

	// calculate stiffness matrix and right part
	grid.assemble_slae(mA, vRightPart);

	// solve SLAE
	size_t nIterationNum = Solvers::CGM(mA, vSolution, vRightPart, 1e-9);

	// calc max error
	double rMaxError = 0.0;
	for (size_t nNode = 0; nNode < grid.nodes_number(); ++nNode)
	{
		rMaxError = std::max(fabs(vSolution[nNode] - Task::exact_solution(grid.coordinates(nNode))), rMaxError);
	}

	std::cout << "Error = " << std::setprecision(10) << rMaxError << std::endl;
	std::cout << "#Iterations = " << nIterationNum << std::endl;

	std::cout << "Solution:\n";
	for (size_t nRow = 0; nRow < nGridParam; ++nRow)
	{
		for (size_t nColumn = 0; nColumn < nGridParam; ++nColumn)
		{
			std::cout << vSolution[nRow * nGridParam + nColumn] << ' ';
		}
		std::cout << std::endl;
	}

	std::cout << "Right part:\n";
	for (size_t nInd = 0; nInd < nGridParam; ++nInd)
	{
		std::cout << vRightPart[nInd] << std::endl;
	}

	return 0;
}


#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "TaskFunctions.h"
#include "FDMApproximation.h"
#include "Solvers.h"

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
	fdm_slau_assembling(grid, mA, vRightPart);

	// solve SLAE
	size_t nIterationNum = PCGM(mA, vSolution, vRightPart, 1e-9);

	// calc max error
	double rMaxError = 0.0;
	for (size_t nNode = 0; nNode < grid.nodes_number(); ++nNode)
	{
		rMaxError = std::max(fabs(vSolution[nNode] - Task::exact_solution(grid.coordinates(nNode))), rMaxError);
	}

	std::cout << "h = " << std::setprecision(5) << grid.m_rH << std::endl;
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

	//int    nI, nJ;
	/*for (nI = 0; nI < mA.n_rows(); nI++)
	{
	  printf("\n");
	  for (nJ = 0; nJ < mA.n_columns(); nJ++)
	  {
		printf("%5.1f ", mA(nI, nJ));
	  }
	}
	printf("\n");
  */
	//sparse_matrix::iterator itA(mA);

	/*for (itA.First1(); !itA.IsDone1(); itA.Next1())
	{
	  printf("\n");
	  for (itA.First2(); !itA.IsDone2(); itA.Next2())
	  {
		printf("%5.1f(%d %d) ", *itA, itA.first_index(), itA.second_index());
	  }
	}*/

	/*for (nI = 0; nI < mA.n_rows(); nI++)
	{
	  printf("\n");
	  for (itA.First2(nI); !itA.IsDone2(); itA.Next2())
	  {
		printf("%5.1f(%d %d) ", *itA, itA.first_index(), itA.second_index());
	  }
	}*/

	return 0;
}


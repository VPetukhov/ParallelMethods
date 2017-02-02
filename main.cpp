#include <stdio.h>

#include "TaskFunctions.h"
#include "FDMApproximation.h"
#include "Solvers.h"

int main()
{
	int nGridParam, nBWidth, nNode;

	// set the partition number of the square side
	nGridParam = 4;
	// set maximal number of nonzero matrix elements in row
	nBWidth = 5;

	// create grid
	fdm_grid grid(nGridParam);

	// allocate SLAE
	// - sparse matrix
	sparse_matrix mA(grid.m_nSize, nBWidth);
	// - solution and right part
	vector vSolution(grid.m_nSize), vRightPart(grid.m_nSize);

	// calculate stiffness matrix and right part
	fdm_slau_assembling(grid, mA, vRightPart);

	// solve SLAE
	PCGM(mA, vSolution, vRightPart, 1e-9);

	// calc max error
	double rMaxError = 0.0;
	for (nNode = 0; nNode < grid.m_nSize; ++nNode)
	{
		rMaxError = std::max(fabs(vSolution[nNode] - exact_solution(grid.m_vCoords[nNode])), rMaxError);
	}

	printf("h     = %1.5f\n", grid.m_rH);
	printf("Error = %1.10f\n", rMaxError);

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

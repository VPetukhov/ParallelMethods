#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <mpi.h>

#include <Utils/FDMGridMPI.h>
#include <Utils/TaskFunctions.h>

#include "Solvers/CGM_MPI.h"
#include "SparseMatrix.h"

int main(int argc, char **argv)
{
	if (int nRc = MPI_Init(&argc, &argv))
	{
		std:: cout << "Error: " << nRc << std::endl;
		MPI_Abort(MPI_COMM_WORLD, nRc);
	}

	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
//	MPI_Comm_size(MPI_COMM_WORLD, &num_threads);

	// create grid
	int nProcVertical, nProcHorizontal;
	double rHeight, rWidth, rDh, rDw;

	size_t nBWidth = 5;
	nProcHorizontal = 5; nProcVertical = 4;
//	nProcHorizontal = 1; nProcVertical = 1;
	rHeight = 1; rWidth = 1; rDh = 0.01, rDw = 0.02;

	fdm_grid_mpi grid(nRank, nProcVertical, nProcHorizontal, rHeight, rWidth, rDh, rDw);
	sparse_matrix mA(grid.nodes_number(), nBWidth);

	vector vSolution(grid.nodes_number()), vRightPart(grid.nodes_number());

	Solvers::cgm_mpi cgm(grid.boundaries());

	grid.assemble_slae(mA, vRightPart);

	size_t nIterationNum = cgm.solve(mA, vSolution, vRightPart, 1e-6);

	if (nRank == 5)
	{
		std::cout << "5 Sizes: " << grid.rows_number() << " " << grid.columns_number() << std::endl;
	}
	if (nRank == 0)
	{
		std::cout << "0 Sizes: " << grid.rows_number() << " " << grid.columns_number() << std::endl;
	}
	// calc max error
	if (nRank == 0)
	{
		double rMaxError = 0.0;
		for (size_t nNode = 0; nNode < grid.nodes_number(); ++nNode)
		{
			rMaxError = std::max(fabs(vSolution[nNode] - Task::exact_solution(grid.coordinates(nNode))), rMaxError);
		}

		for (size_t nRow = 0; nRow < grid.rows_number(); ++nRow)
		{
			for (size_t nColumn = 0; nColumn < grid.columns_number(); ++nColumn)
			{
				size_t nNode = nRow * grid.columns_number() + nColumn;
//				std::cout << std::setprecision(3) << std::round(fabs(vSolution[nNode] - Task::exact_solution(grid.coordinates(nNode))) * 1000) / 1000 << '\t';
				std::cout << std::setprecision(3) << vSolution[nNode] << '\t';
//				std::cout << std::setprecision(3) << Task::exact_solution(grid.coordinates(nNode)) << '\t';
			}
			std::cout << std::endl;
		}

		std::cout << "Error = " << std::setprecision(10) << rMaxError << std::endl;
		std::cout << "#Iterations = " << nIterationNum << std::endl;
	}

//	std::cout << "Solution:\n";
//	for (size_t nRow = 0; nRow < nGridParam; ++nRow)
//	{
//		for (size_t nColumn = 0; nColumn < nGridParam; ++nColumn)
//		{
//			std::cout << vSolution[nRow * nGridParam + nColumn] << ' ';
//		}
//		std::cout << std::endl;
//	}
//
//	std::cout << "Right part:\n";
//	for (size_t nInd = 0; nInd < nGridParam; ++nInd)
//	{
//		std::cout << vRightPart[nInd] << std::endl;
//	}

	MPI_Finalize();
	return 0;
}


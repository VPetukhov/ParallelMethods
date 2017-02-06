#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <mpi.h>

#include <Utils/FDMGridMPI.h>
#include <Utils/TaskFunctions.h>

#include "Solvers/CGM_MPI.h"
#include "SparseMatrix.h"

void get_proc_number(int &nProcHorizontal, int &nProcVertical)
{
	int nThreadsNum;
	MPI_Comm_size(MPI_COMM_WORLD, &nThreadsNum);

	int sqr_nthreads = (int)std::sqrt(nThreadsNum);
	while (sqr_nthreads > 0)
	{
		if (nThreadsNum % sqr_nthreads == 0)
		{
			nProcVertical = nThreadsNum / sqr_nthreads;
			nProcHorizontal = nThreadsNum / nProcVertical;
			sqr_nthreads = 0;
		}
		else sqr_nthreads--;
	}
}

template <typename T>
std::ostream& print_vector(const T &vArr, int nIndDisplace = 0, int nMaxInd = -1, const std::string &sPrefix = "",
                           bool bOutIndex=false, const std::string &sDelim = "\t")
{
	if (nMaxInd < 0)
	{
		nMaxInd = vArr.size() - nIndDisplace;
	}
	std::ostringstream str;
	if (sPrefix != "")
	{
		str << sPrefix << sDelim;
	}

	if (bOutIndex)
	{
		for (int i = 0; i < nMaxInd; ++i)
		{
			str << nIndDisplace + i << sDelim;
		}
	}
	else
	{
		for (int i = 0; i < nMaxInd; ++i)
		{
			str << vArr[nIndDisplace + i] << sDelim;
		}
	}

	return std::cout << str.str() << std::endl;
}

std::vector<double> gather_solution(vector &vSolutionPart)
{
	int nProcessNumber;
	int nSolutionSize = vSolutionPart.size();
	MPI_Comm_size(MPI_COMM_WORLD, &nProcessNumber);

	std::vector<int> vSizes(nProcessNumber), vDispls;
	MPI_Allgather(&nSolutionSize, 1, MPI_INT, vSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
	int nFullSolutionSize = 0;
	for (int nSize : vSizes)
	{
		vDispls.push_back(nFullSolutionSize);
		nFullSolutionSize += nSize;
	}

	std::vector<double> vSolution(nFullSolutionSize);
	MPI_Gatherv(&vSolutionPart[0], nSolutionSize, MPI_DOUBLE, vSolution.data(), vSizes.data(), vDispls.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return vSolution;
}

void print_solution(const std::vector<double> &vSolution, double rHeight, double rWidth, double rDh, double rDw)
{
	int nColumnsNumber = (int) (rWidth / rDw + 0.5);
	int nRow = 0;
	std::cout << "Solution:\n";
	for (double rY = 0; rY < rHeight; rY += rDh, ++nRow)
	{
		print_vector(vSolution, nRow * nColumnsNumber + 1, nColumnsNumber - 2);
	}

	std::cout << "\nError:\n";
	nRow = 1;
	double rMaxError = 0.0;
	for (double rY = rDh; rY < rHeight + rDh; rY += rDh, ++nRow)
	{
		int nColumn = 1;
		for (double rX = rDw; rX < rWidth - rDw; rX += rDw, ++nColumn)
		{
			double rCurError = std::abs(vSolution[nRow * nColumnsNumber + nColumn] - Task::exact_solution(fdm_grid_mpi::coord(rX, rY)));
			rMaxError = std::max(rCurError, rMaxError);
			std::cout << rCurError << "\t";
		}
		std::cout << "\n";
	}

	std::cout << "Max error: " << rMaxError << std::endl;
}

int main(int argc, char **argv)
{
	if (int nRc = MPI_Init(&argc, &argv))
	{
		std:: cout << "Error: " << nRc << std::endl;
		MPI_Abort(MPI_COMM_WORLD, nRc);
	}

	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	double rHeight, rWidth, rDh, rDw;
	int nProcVertical, nProcHorizontal;
	get_proc_number(nProcHorizontal, nProcVertical);

	size_t nBWidth = 5;
	rHeight = 1; rWidth = 1; rDh = 0.1, rDw = 0.1;

	fdm_grid_mpi grid(nRank, nProcVertical, nProcHorizontal, rHeight, rWidth, rDh, rDw);
	sparse_matrix mA(grid.nodes_number(), nBWidth);

	vector vSolution(grid.nodes_number()), vRightPart(grid.nodes_number());

	Solvers::cgm_mpi cgm(grid.boundaries(), grid.neighbours(), grid.bound_flags());

	grid.assemble_slae(mA, vRightPart);

	size_t nIterationNum = cgm.solve(mA, vSolution, vRightPart, 1e-6, (int) (rHeight / rDh * rWidth / rDw) * 5);

	// calc max error
	std::vector<double> vSolutionFull = gather_solution(vSolution);
	if (nRank == 0)
	{
		print_solution(vSolutionFull, rHeight, rWidth, rDh, rDw);
	}

	MPI_Finalize();
	return 0;
}


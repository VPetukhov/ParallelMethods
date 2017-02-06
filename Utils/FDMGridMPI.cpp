#include "FDMGridMPI.h"

#include "SparseMatrix.h"
#include "TaskFunctions.h"

#include <iostream>
#include <stdexcept>
#include <mpi.h>

fdm_grid_mpi::fdm_grid_mpi(int nRank, int nProcVertical, int nProcHorizontal, double rHeight, double rWidth,
                           double rDh, double rDw)
	: m_vBoundaries(SIDES_NUMBER)
{
	int nProcRow = nRank / nProcHorizontal;
	int nProcColumn = nRank % nProcHorizontal;
	double rStartX = rWidth / nProcHorizontal * nProcColumn;
	double rEndX = rWidth / nProcHorizontal * (nProcColumn + 1);
	double rStartY = rHeight / nProcVertical * nProcRow;
	double rEndY = rHeight / nProcVertical * (nProcRow + 1);

	if (nProcColumn > 0)
	{
		rStartX -= rDw;
	}

	if (nProcRow > 0)
	{
		rStartY -= rDh;
	}

	std::vector<int> vBoundaryProcesses(SIDES_NUMBER);
	vBoundaryProcesses[TOP] = nProcRow == (nProcVertical - 1) ? NOT_PROCESS : nProcHorizontal * (nProcRow + 1) + nProcColumn;
	vBoundaryProcesses[BOTTOM] = nProcRow == 0 ? NOT_PROCESS : nProcHorizontal * (nProcRow - 1) + nProcColumn;
	vBoundaryProcesses[LEFT] = nProcColumn == 0 ? NOT_PROCESS : nRank - 1;
	vBoundaryProcesses[RIGHT] = nProcColumn == (nProcHorizontal - 1) ? NOT_PROCESS : nRank + 1;

//	std::cout << nRank << ": (" << rStartX << ", " << rStartY << "); (" << rEndX << ", " << rEndY << "); ("
//	          << rDw << ", " << rDh << "); " << std::endl;
	this->init(vBoundaryProcesses, rStartX, rStartY, rEndX, rEndY, rDw, rDh);
}

void fdm_grid_mpi::init(std::vector<int> vBoundaryProcesses, double rStartX, double rStartY, double rEndX, double rEndY,
                        double rDw, double rDh)
{
	this->m_nRowsNum = (size_t) ((rEndY - rStartY) / rDh + 0.5) + 1;
	this->m_nColumnsNum = (size_t) ((rEndX - rStartX) / rDw + 0.5) + 1;

	this->m_nNodesNumber = this->m_nRowsNum * this->m_nColumnsNum;

	this->m_rDS = rDw * rDh;

	for (int nSideInd = 0; nSideInd < SIDES_NUMBER; ++nSideInd)
	{
		this->m_vBoundaries[nSideInd].nProcessRank = vBoundaryProcesses[nSideInd];
	}

	this->m_vIsBound.resize(this->m_nNodesNumber, false);
	this->m_vIsConstantBound.resize(this->m_nNodesNumber, false);

	for (size_t nRow = 0; nRow < this->m_nRowsNum; ++nRow)
	{
		// set coordinates
		double rYCoord = rStartY + nRow * rDh;
		for (size_t nColumn = 0; nColumn < this->m_nColumnsNum; ++nColumn)
		{
			this->m_vCoords.push_back(coord(rStartX + nColumn * rDw,  rYCoord));
			this->m_vNeighbours.push_back(this->fill_neighbours(nRow, nColumn, this->m_nRowsNum, this->m_nColumnsNum));
		}
	}

	// Boundaries
	size_t nLastRowInd = (this->m_nRowsNum - 1) * this->m_nColumnsNum;
	for (size_t nColumn = 1; nColumn < this->m_nColumnsNum - 1; ++nColumn)
	{
		this->set_boundary(nColumn, nColumn + this->m_nColumnsNum, BOTTOM);
		this->set_boundary(nLastRowInd + nColumn, nLastRowInd + nColumn - this->m_nColumnsNum, TOP);
	}

	for (size_t nRow = 1; nRow < this->m_nRowsNum - 1; ++nRow)
	{
		size_t nLineInd = nRow * this->m_nColumnsNum;
		this->set_boundary(nLineInd, nLineInd + 1, LEFT);
		this->set_boundary(nLineInd + this->m_nColumnsNum - 1, nLineInd + this->m_nColumnsNum - 2, RIGHT);
	}

	this->m_vRedundantNodes.push_back(0);
	this->m_vRedundantNodes.push_back(this->m_nColumnsNum - 1);
	this->m_vRedundantNodes.push_back(nLastRowInd);
	this->m_vRedundantNodes.push_back(nLastRowInd + this->m_nColumnsNum - 1);

//	int nRank;
//	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
//
//	for (int i = 0; i < side_index::SIDES_NUMBER; ++i)
//	{
//		std::cout << nRank << " (" << this->m_vBoundaries[i].nProcessRank << "): ";
//		for (int receive : this->m_vBoundaries[i].vBoundariesToReceive)
//		{
//			std::cout << "(" << this->m_vCoords[receive].rX << ", " << this->m_vCoords[receive].rY << ") ";
//		}
//		std::cout << std::endl;
//	}
}

void fdm_grid_mpi::set_boundary(size_t nRecieveIndex, size_t nSendIndex, side_index nSideIndex)
{
	this->m_vIsBound[nRecieveIndex] = true;

	if (this->m_vBoundaries[nSideIndex].nProcessRank == NOT_PROCESS)
	{
		this->m_vIsConstantBound[nRecieveIndex] = true;
	}
	else
	{
		this->m_vBoundaries[nSideIndex].vBoundariesToSend.push_back(nSendIndex);
		this->m_vBoundaries[nSideIndex].vBoundariesToReceive.push_back(nRecieveIndex);
	}
}

std::vector<size_t> fdm_grid_mpi::fill_neighbours(size_t nRow, size_t nColumn, size_t nRowsNum, size_t nColumnsNum)
{
	std::vector<size_t> result(SIDES_NUMBER);
	result[BOTTOM] = nRow == (nRowsNum - 1) ? EMPTY_NEIGHBOUR : (nRow + 1) * nColumnsNum + nColumn;
	result[TOP] = nRow == 0 ? EMPTY_NEIGHBOUR : (nRow - 1) * nColumnsNum + nColumn;
	result[LEFT] = nColumn == 0 ? EMPTY_NEIGHBOUR : nRow * nColumnsNum + nColumn - 1;
	result[RIGHT] = nColumn == (nColumnsNum - 1) ? EMPTY_NEIGHBOUR : nRow * nColumnsNum + nColumn + 1;

	if (nColumn == nColumnsNum || nColumn == 0)
	{
		if (nRow == nRowsNum - 1)
		{
			result[TOP] = EMPTY_NEIGHBOUR;
		}
		else if (nRow == 1)
		{
			result[BOTTOM] = EMPTY_NEIGHBOUR;
		}
	}

	if (nRow == nRowsNum || nRow == 0)
	{
		if (nColumn == nColumnsNum - 1)
		{
			result[RIGHT] = EMPTY_NEIGHBOUR;
		}
		else if (nColumn == 1)
		{
			result[LEFT] = EMPTY_NEIGHBOUR;
		}
	}

	return result;
}

void fdm_grid_mpi::assemble_slae(sparse_matrix &mSM, vector &vRightPart) const
{
	auto redundantNodeIter = this->m_vRedundantNodes.begin();
	for (size_t nNode = 0; nNode < this->nodes_number(); ++nNode)
	{
		if (redundantNodeIter != this->m_vRedundantNodes.end() && *redundantNodeIter == nNode)
		{
			++redundantNodeIter;
			continue;
		}

		if (this->m_vIsBound[nNode])
		{
			// write boundary condition
			mSM.at(nNode, nNode) = 1.0;
			vRightPart[nNode] = this->m_vIsConstantBound[nNode] ? Task::bound_condition(this->m_vCoords[nNode]) : 0;
		}
		else
		{
			// write finite difference equation
			vRightPart[nNode] = this->m_rDS * Task::right_part(this->m_vCoords[nNode]);
			mSM.at(nNode, nNode) = 4.0;

			for (size_t neighbourInd : this->m_vNeighbours[nNode])
			{
				if (neighbourInd == fdm_grid_mpi::EMPTY_NEIGHBOUR)
					throw std::runtime_error("Unexpected index for neighbour to " + std::to_string(nNode));

				if (neighbourInd > this->m_vCoords.size())
					throw std::runtime_error("Index for neighbour to " + std::to_string(nNode) + " is too large (" +
							                         std::to_string(neighbourInd) + " of " + std::to_string(this->m_vCoords.size()) + ")!");

				if (!this->m_vIsBound[neighbourInd])
				{
					mSM.at(nNode, neighbourInd) = -1.0;
				}
				else if (this->m_vIsConstantBound[neighbourInd])
				{
					vRightPart[nNode] -= -1.0 * Task::bound_condition(this->m_vCoords[neighbourInd]);
				}
			}
		}
	}
}

const std::vector<fdm_grid_mpi::boundaries_list> &fdm_grid_mpi::boundaries() const
{
	return this->m_vBoundaries;
}

size_t fdm_grid_mpi::nodes_number() const
{
	return this->m_nNodesNumber;
}

const fdm_grid_mpi::coord &fdm_grid_mpi::coordinates(size_t nIndex) const
{
	return this->m_vCoords.at(nIndex);
}

size_t fdm_grid_mpi::rows_number() const
{
	return this->m_nRowsNum;
}

size_t fdm_grid_mpi::columns_number() const
{
	return this->m_nColumnsNum;
}

const std::vector<std::vector<size_t>> &fdm_grid_mpi::neighbours() const
{
	return this->m_vNeighbours;
}

const std::vector<bool> &fdm_grid_mpi::bound_flags() const
{
	return this->m_vIsBound;
}

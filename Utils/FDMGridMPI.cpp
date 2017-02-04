#include "FDMGridMPI.h"

#include "SparseMatrix.h"
#include "TaskFunctions.h"

#include <iostream>
#include <stdexcept>

fdm_grid_mpi::fdm_grid_mpi(int nRank, int nProcVertical, int nProcHorizontal, double rHeight, double rWidth,
                           double rDh, double rDw)
	: m_vBoundaries(SIDES_NUMBER)
{
	int nRow = nRank / nProcHorizontal;
	int nColumn = nRank % nProcHorizontal;
	double rStartX = rWidth / nProcHorizontal * nColumn;
	double rEndX = rWidth / nProcHorizontal * (nColumn + 1);
	double rStartY = rHeight / nProcVertical * nRow;
	double rEndY = rHeight / nProcVertical * (nRow + 1);

	std::vector<int> vBoundaryProcesses(SIDES_NUMBER);
	vBoundaryProcesses[LEFT] = nColumn == 0 ? NOT_PROCESS : nRank - 1;
	vBoundaryProcesses[RIGHT] = nColumn == nProcHorizontal ? NOT_PROCESS : nRank + 1;
	vBoundaryProcesses[BOTTOM] = nRow == nProcVertical ? NOT_PROCESS : nProcHorizontal * (nRow - 1) + nColumn;
	vBoundaryProcesses[TOP] = nRow == 0 ? NOT_PROCESS : nProcHorizontal * (nRow + 1) + nColumn;

	this->init((size_t) ((rEndX - rStartX) / rDw), (size_t) ((rEndY - rStartY) / rDh), vBoundaryProcesses, rStartX, rStartY, rDw, rDh);
}

void fdm_grid_mpi::init(size_t nRowsNum, size_t nColumnsNum, std::vector<int> vBoundaryProcesses, double rStartX,
                        double rStartY, double rDw, double rDh)
{
	nRowsNum += 2;
	nColumnsNum += 2;

	this->m_nNodesNumber = (nRowsNum + 1) * (nColumnsNum + 1);
	this->m_rDS = rDw * rDh;

	for (int nSideInd = 0; nSideInd < SIDES_NUMBER; ++nSideInd)
	{
		this->m_vBoundaries[nSideInd].process_rank = vBoundaryProcesses[nSideInd];
	}

	this->m_vIsBound.resize(this->m_nNodesNumber, false);
	this->m_vIsConstantBound.resize(this->m_nNodesNumber, false);

	for (size_t nRow = 0; nRow < nRowsNum; ++nRow)
	{
		// set coordinates
		double rYCoord = rStartY + nRow * rDh;
		for (size_t nColumn = 0; nColumn < nColumnsNum; ++nColumn)
		{
			this->m_vCoords.push_back(coord(rStartX + nColumn * rDw,  rYCoord));
			this->m_vNeighbours.push_back(this->fill_neighbours(nRow, nColumn, nRowsNum, nColumnsNum));
		}

		// - horizontal boundary
		size_t nLineInd = nRow * nColumnsNum;
		this->set_boundary(nLineInd, nLineInd + 1, LEFT);
		this->set_boundary(nLineInd + nColumnsNum - 1, nLineInd + nColumnsNum - 2, RIGHT);
	}

	// - vertical boundary
	size_t nLastRowInd = (nRowsNum - 1) * nRowsNum;
	for (size_t nColumn = 0; nColumn < nColumnsNum; ++nColumn)
	{
		this->set_boundary(nColumn, nColumn + nColumnsNum, TOP);
		this->set_boundary(nLastRowInd + nColumn, nLastRowInd + nColumn - nColumnsNum, BOTTOM);
	}
}

void fdm_grid_mpi::set_boundary(size_t nRecieveIndex, size_t nSendIndex, side_index nSideIndex)
{
	this->m_vIsBound[nRecieveIndex] = true;

	if (this->m_vBoundaries[nSideIndex].process_rank == NOT_PROCESS)
	{
		this->m_vIsConstantBound[nRecieveIndex] = true;
	}
	else
	{
		this->m_vBoundaries[nSideIndex].boundaries_to_send.push_back(nSendIndex);
		this->m_vBoundaries[nSideIndex].boundaries_to_receive.push_back(nRecieveIndex);
	}
}

std::vector<size_t> fdm_grid_mpi::fill_neighbours(size_t nRow, size_t nColumn, size_t nRowsNum, size_t nColumnsNum)
{
	std::vector<size_t> result(SIDES_NUMBER);
	result[BOTTOM] = nRow == nRowsNum ? EMPTY_NEIGHBOUR : (nRow + 1) * nColumnsNum + nColumn;
	result[TOP] = nRow == 0 ? EMPTY_NEIGHBOUR : (nRow - 1) * nColumnsNum + nColumn;
	result[LEFT] = nColumn == 0 ? EMPTY_NEIGHBOUR : nRow * nColumnsNum + nColumn - 1;
	result[RIGHT] = nColumn == nColumnsNum ? EMPTY_NEIGHBOUR : nRow * nColumnsNum + nColumn + 1;
	return result;
}

void fdm_grid_mpi::assemble_slae(sparse_matrix &mSM, vector &vRightPart) const
{
	for (size_t nNode = 0; nNode < this->nodes_number(); ++nNode)
	{
		if (this->m_vIsBound[nNode])
		{
			// write boundary condition
			mSM.at(nNode, nNode) = 1.0;
			vRightPart[nNode] = this->m_vIsConstantBound[nNode] ? Task::bound_condition(this->coordinates(nNode)) : 0;
		}
		else
		{
			// write finite difference equation
			vRightPart[nNode] = this->m_rDS * Task::right_part(this->coordinates(nNode));
			mSM.at(nNode, nNode) = 4.0;

			for (size_t neighbourInd : this->m_vNeighbours[nNode])
			{
				if (neighbourInd == fdm_grid_mpi::EMPTY_NEIGHBOUR)
					throw std::runtime_error("Unexpected index for neighbour to " + std::to_string(nNode));

				if (!this->m_vIsBound[neighbourInd])
				{
					mSM.at(nNode, neighbourInd) = -1.0;
				}
				else if (this->m_vIsConstantBound[neighbourInd])
				{
					vRightPart[nNode] -= -1.0 * Task::bound_condition(this->coordinates(neighbourInd));
				}
			}
		}
	}
}

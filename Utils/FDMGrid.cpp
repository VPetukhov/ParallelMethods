#include "FDMGrid.h"

#include "SparseMatrix.h"
#include "TaskFunctions.h"

#include <iostream>
#include <stdexcept>

fdm_grid::fdm_grid(size_t nParam)
	: m_rH(1.0 / nParam)
	, m_nNodesNumber((nParam + 1) * (nParam + 1))
{
	size_t nSideSize = nParam + 1;
	this->m_vBoundFlag.resize(this->m_nNodesNumber, false);

	for (size_t nRow = 0; nRow < nSideSize; ++nRow)
	{
		// set coordinates
		size_t nXLine = nRow * nSideSize;
		double rXCoord = nRow * this->m_rH;
		for (size_t nColumn = 0; nColumn < nSideSize; ++nColumn)
		{
			this->m_vCoords.push_back(coord(rXCoord, nColumn * this->m_rH));
			this->m_vNeighbours.push_back(this->fill_neighbours(nRow, nColumn, nSideSize, nSideSize));
		}

		// mark boundary nodes
		// - horizontal boundary
		this->m_vBoundFlag[nXLine] = true;
		this->m_vBoundFlag[nXLine + nSideSize - 1] = true;
		// - vertical boundary
		this->m_vBoundFlag[0 + nRow] = true;
		this->m_vBoundFlag[(nSideSize - 1) * nSideSize + nRow] = true;
	}
}

std::vector<size_t> fdm_grid::fill_neighbours(size_t nRow, size_t nColumn, size_t nRowsNum, size_t nColumnsNum)
{
	std::vector<size_t> result;
	result.push_back(nRow == nRowsNum ? EMPTY_NEIGHBOUR : (nRow + 1) * nColumnsNum + nColumn);
	result.push_back(nRow == 0 ? EMPTY_NEIGHBOUR : (nRow - 1) * nColumnsNum + nColumn);
	result.push_back(nColumn == 0 ? EMPTY_NEIGHBOUR : nRow * nColumnsNum + nColumn - 1);
	result.push_back(nColumn == nColumnsNum ? EMPTY_NEIGHBOUR : nRow * nColumnsNum + nColumn + 1);
	return result;
}

void fdm_grid::assemble_slae(sparse_matrix &mSM, vector &vRightPart) const
{
	double rH2 = this->m_rH * this->m_rH;

	// assemble SLAE
	for (size_t nNode = 0; nNode < this->nodes_number(); ++nNode)
	{
		if (this->node_is_boundary(nNode))
		{
			// write boundary condition
			vRightPart[nNode] = Task::bound_condition(this->coordinates(nNode));
			mSM.at(nNode, nNode) = 1.0;
		}
		else
		{
			// write finite difference equation
			vRightPart[nNode] = rH2 * Task::right_part(this->coordinates(nNode));
			mSM.at(nNode, nNode) = 4.0;

			for (size_t neighbourInd : this->m_vNeighbours[nNode])
			{
				if (neighbourInd == fdm_grid::EMPTY_NEIGHBOUR)
					throw std::runtime_error("Unexpected index for neighbour to " + std::to_string(nNode));

				if (this->node_is_boundary(neighbourInd))
				{
					// take into account boundary condition for SLAE symmetrization
					vRightPart[nNode] -= -1.0 * Task::bound_condition(this->coordinates(neighbourInd));
				}
				else
				{
					mSM.at(nNode, neighbourInd) = -1.0;
				}
			}
		}
	}
}

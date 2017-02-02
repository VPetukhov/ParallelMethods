#include <iostream>
#include <stdexcept>

#include "FDMGrid.h"

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


const std::vector<size_t>& fdm_grid::neighbours(size_t nNode) const
{
	return this->m_vNeighbours.at(nNode);
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

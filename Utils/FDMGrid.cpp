#include <iostream>
#include <stdexcept>

#include "FDMGrid.h"

fdm_grid::fdm_grid(size_t nParam)
	: m_rH(1.0 / nParam)
	, m_nN(nParam + 1)
	, m_nSize((nParam + 1) * (nParam + 1))
{
	this->m_vCoords.resize(this->m_nSize);
	this->m_vBoundFlag.resize(this->m_nSize, false);


	for (size_t nI = 0; nI < this->m_nN; ++nI)
	{
		// set coordinates
		size_t nXLine = nI * this->m_nN;
		double rXCoord = nI * this->m_rH;
		for (size_t nJ = 0; nJ < this->m_nN; ++nJ)
		{
			this->m_vCoords[nXLine + nJ].rX = rXCoord;
			this->m_vCoords[nXLine + nJ].rY = nJ * this->m_rH;
		}

		// mark boundary nodes
		// - horizontal boundary
		this->m_vBoundFlag[nXLine] = true;
		this->m_vBoundFlag[nXLine + this->m_nN - 1] = true;
		// - vertical boundary
		this->m_vBoundFlag[0 + nI] = true;
		this->m_vBoundFlag[(this->m_nN - 1) * this->m_nN + nI] = true;
	}
}


void fdm_grid::fill_neighbours(size_t nNode, std::vector<size_t> &vNeighbours) const
{
	if (vNeighbours.size() != 4)
		throw std::runtime_error("invalid neighbours vector size");

	vNeighbours[0] = nNode - 1;
	vNeighbours[1] = nNode + 1;
	vNeighbours[2] = nNode - this->m_nN;
	vNeighbours[3] = nNode + this->m_nN;
}


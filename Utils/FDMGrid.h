#pragma once

#include <vector>
#include <cstddef>

struct coord
{
	double rX, rY;
};

class fdm_grid
{
private:
	std::vector<coord> m_vCoords;
	std::vector<bool> m_vBoundFlag;

public:
	const size_t m_nN, m_nSize;
	const double m_rH;

public:
	fdm_grid(size_t nParam);

	~fdm_grid()
	{}

	bool node_is_boundary(size_t nNode) const
	{
		return this->m_vBoundFlag.at(nNode);
	}

	const coord& coordinates(size_t index) const
	{
		return this->m_vCoords.at(index);
	}

	void fill_neighbours(size_t nNode, std::vector<size_t> &vNeighbours) const;
};

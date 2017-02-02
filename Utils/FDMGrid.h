#pragma once

#include <cstddef>
#include <vector>
#include <limits>

class fdm_grid
{
public:
	struct coord
	{
		double rX, rY;

		coord(double rX = 0, double rY = 0)
			: rX(rX), rY(rY)
		{}
	};

private:
	std::vector<coord> m_vCoords;
	std::vector<std::vector<size_t>> m_vNeighbours;
	std::vector<bool> m_vBoundFlag;

	const size_t m_nNodesNumber;
public:
	static const size_t EMPTY_NEIGHBOUR = std::numeric_limits<size_t>::max();
	const double m_rH;

private:
	std::vector<size_t> fill_neighbours(size_t nRow, size_t nColumn, size_t nRowsNum, size_t nColumnsNum);

public:
	fdm_grid(size_t nParam);

	~fdm_grid()
	{}

	bool node_is_boundary(size_t nNode) const
	{
		return this->m_vBoundFlag.at(nNode);
	}

	const coord& coordinates(size_t nIndex) const
	{
		return this->m_vCoords.at(nIndex);
	}

	size_t nodes_number() const
	{
		return this->m_nNodesNumber;
	}

	const std::vector<size_t>& neighbours(size_t nNode) const;
};

#pragma once

#include <cstddef>
#include <vector>
#include <limits>
#include "FDMGrid.h"

class sparse_matrix;
class vector;

class fdm_grid_mpi
{
public:
	typedef fdm_grid::coord coord;

	struct boundaries_list
	{
		std::vector<int> boundaries_to_send;
		std::vector<int> boundaries_to_receive;
		int process_rank;
	};

	enum side_index : int
	{
		TOP,
		BOTTOM,
		LEFT,
		RIGHT,
		SIDES_NUMBER
	};

	static const int NOT_PROCESS = -1;

private:
	std::vector<coord> m_vCoords;
	std::vector<std::vector<size_t>> m_vNeighbours;
	std::vector<bool> m_vIsBound, m_vIsConstantBound;
	std::vector<boundaries_list> m_vBoundaries;

	size_t m_nNodesNumber;
	double m_rDS;

	static const size_t EMPTY_NEIGHBOUR = std::numeric_limits<size_t>::max();

private:
	std::vector<size_t> fill_neighbours(size_t nRow, size_t nColumn, size_t nRowsNum, size_t nColumnsNum);

	void init(size_t nRowsNum, size_t nColumnsNum, std::vector<int> vBoundaryProcesses, double rStartX,
	          double rStartY, double rDw, double rDh);

public:
	fdm_grid_mpi(int nRank, int nProcVertical, int nProcHorizontal, double rHeight, double rWidth,
		             double rDh, double rDw);

	~fdm_grid_mpi()
	{}

	const coord& coordinates(size_t nIndex) const
	{
		return this->m_vCoords.at(nIndex);
	}

	size_t nodes_number() const
	{
		return this->m_nNodesNumber;
	}

	void assemble_slae(sparse_matrix &mSM, vector &vRightPart) const;

	void set_boundary(size_t nRecieveIndex, size_t nSendIndex, side_index nSideIndex);
};

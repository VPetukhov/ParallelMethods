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
		std::vector<int> vBoundariesToSend;
		std::vector<int> vBoundariesToReceive;
		int nProcessRank;
	};

	enum side_index : int
	{
		TOP = 0,
		BOTTOM,
		LEFT,
		RIGHT,
		SIDES_NUMBER
	};

	static const int NOT_PROCESS = -1;
	static const size_t EMPTY_NEIGHBOUR = std::numeric_limits<size_t>::max();

private:
	std::vector<coord> m_vCoords;
	std::vector<std::vector<size_t>> m_vNeighbours;
	std::vector<bool> m_vIsBound, m_vIsConstantBound;
	std::vector<int> m_vRedundantNodes;
	std::vector<boundaries_list> m_vBoundaries;

	size_t m_nNodesNumber, m_nRowsNum, m_nColumnsNum;
	double m_rDS;

private:
	std::vector<size_t> fill_neighbours(size_t nRow, size_t nColumn, size_t nRowsNum, size_t nColumnsNum);

	void init(std::vector<int> vBoundaryProcesses, double rStartX, double rStartY, double rEndX, double rEndY,
		          double rDw, double rDh);

	void set_boundary(size_t nRecieveIndex, size_t nSendIndex, side_index nSideIndex);

public:
	fdm_grid_mpi(int nRank, int nProcVertical, int nProcHorizontal, double rHeight, double rWidth,
		             double rDh, double rDw);

	const coord& coordinates(size_t nIndex) const;
	size_t nodes_number() const;
	size_t rows_number() const;
	size_t columns_number() const;

	const std::vector<boundaries_list> &boundaries() const;
	const std::vector<std::vector<size_t>> &neighbours() const;
	const std::vector<bool> &bound_flags() const;
	void assemble_slae(sparse_matrix &mSM, vector &vRightPart) const;
};

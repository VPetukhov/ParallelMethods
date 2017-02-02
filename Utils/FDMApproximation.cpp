#include "TaskFunctions.h"

#include "FDMApproximation.h"

void fdm_slau_assembling(const fdm_grid &grid, sparse_matrix &mSM, vector &vRightPart)
{
	std::vector<size_t> vNeighbours(4);
	double rH2 = grid.m_rH * grid.m_rH;

	// assemble SLAE
	for (size_t nNode = 0; nNode < grid.m_nSize; ++nNode)
	{
		if (grid.node_is_boundary(nNode))
		{
			// write boundary condition
			vRightPart[nNode] = Task::bound_condition(grid.coordinates(nNode));
			mSM.at(nNode, nNode) = 1.0;
		}
		else
		{
			// write finite difference equation
			vRightPart[nNode] = rH2 * Task::right_part(grid.coordinates(nNode));
			mSM.at(nNode, nNode) = 4.0;

			grid.fill_neighbours(nNode, vNeighbours);
			for (size_t neighbourInd : vNeighbours)
			{
				if (grid.node_is_boundary(neighbourInd))
				{
					// take into account boundary condition for SLAE symmetrization
					vRightPart[nNode] -= -1.0 * Task::bound_condition(grid.coordinates(neighbourInd));
				}
				else
				{
					mSM.at(nNode, neighbourInd) = -1.0;
				}
			}
		}
	}
} // end of fdm_slau_assembling


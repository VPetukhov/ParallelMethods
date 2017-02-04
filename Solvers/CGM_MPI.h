#pragma once

#include <cstdlib>
#include "Utils/FDMGridMPI.h"
#include <mpi.h>

class sparse_matrix;
class vector;

namespace Solvers
{
	class cgm_mpi
	{
	private:
		typedef fdm_grid_mpi::side_index side_index;
		struct boundaries_list
		{
			int nProcessRank;
			MPI_Datatype send_t;
			MPI_Datatype receive_t;
		};

	private:
		std::vector<boundaries_list> vBoundaryProcs;

	private:
		std::vector<MPI_Request> update_boundaries(vector &vSolution) const;
		void create_datatype(std::vector<int> indexes, MPI_Datatype *datatype) const;
		side_index get_opposite_side(side_index side) const;

	public:
		cgm_mpi(const std::vector<fdm_grid_mpi::boundaries_list> &vBoundaryProcs);
		size_t solve(const sparse_matrix &mA, vector &vSol, const vector &vRightPart, double rEps);
	};
}

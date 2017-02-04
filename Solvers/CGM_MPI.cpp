#include <mpi.h>
#include "CGM_MPI.h"

#include "SparseMatrix.h"

namespace Solvers
{
	cgm_mpi::cgm_mpi(const std::vector<fdm_grid_mpi::boundaries_list> &vBoundaryProcs)
		: vBoundaryProcs(vBoundaryProcs.size())
	{
		for (int i = 0; i < side_index::SIDES_NUMBER; ++i)
		{
			auto const &boundInfo = vBoundaryProcs[i];
			this->vBoundaryProcs[i].nProcessRank = boundInfo.process_rank;
			this->vBoundaryProcs[i].receive_t = this->create_datatype(boundInfo.boundaries_to_receive);
			this->vBoundaryProcs[i].send_t = this->create_datatype(boundInfo.boundaries_to_send);
		}
	}

	size_t cgm_mpi::solve(const sparse_matrix &mA, vector &vSol, const vector &vRightPart, double rEps)
	{
		size_t  nSize = vRightPart.size();
		vSol = 0.0;

		vector vGradient = vRightPart - mA.vector_multiply(vSol);
		vector vResiduals = vGradient;
		double rBRR1 = vGradient * vGradient;

		size_t nIterationNum;
		std::vector<MPI_Request> vSendRequests(this->vBoundaryProcs.size(), MPI_REQUEST_NULL);
		std::vector<MPI_Status> vSendStatuses(vSendRequests.size());

		for (nIterationNum = 0; (nIterationNum < 2 * nSize) && (vGradient.norm_inf() > rEps); ++nIterationNum)
		{
			vector vTemp = mA.vector_multiply(vResiduals);
			double rAlpha = rBRR1 / (vTemp * vResiduals);

			MPI_Waitall((int)vSendRequests.size(), vSendRequests.data(), vSendStatuses.data());
			vSol.add_scale_vector(vResiduals, rAlpha);

			vSendRequests = this->update_boundaries(vSol);
			vGradient = vRightPart - mA.vector_multiply(vSol); // TODO: Optimize (try to return simple implementation)

			double rBRR2  = vGradient * vGradient;
			double rBetta = rBRR2 / rBRR1;
			rBRR1  = rBRR2;

			vTemp = vGradient;
			vTemp.add_scale_vector(vResiduals, rBetta);
			vResiduals = vTemp;
		}

		return nIterationNum;
	}

	std::vector<MPI_Request> cgm_mpi::update_boundaries(vector &vSolution) const
	{
		std::vector<MPI_Request> vSendRequests(this->vBoundaryProcs.size());
		std::vector<MPI_Request> vReceiveRequests(this->vBoundaryProcs.size());

		for (int i = 0; i < this->vBoundaryProcs.size(); ++i)
		{
			auto const &proc = this->vBoundaryProcs[i];
			MPI_Irecv(&vSolution[0], 1, proc.receive_t, proc.nProcessRank, i, MPI_COMM_WORLD, &vReceiveRequests[i]);
		}

		for (int i = 0; i < this->vBoundaryProcs.size(); ++i)
		{
			auto const &proc = this->vBoundaryProcs[i];
			MPI_Isend(&vSolution[0], 1, proc.send_t, proc.nProcessRank, this->get_opposite_side((side_index)i), MPI_COMM_WORLD, &vSendRequests[i]);
		}

		std::vector<MPI_Status> vReceiveStatuses(vReceiveRequests.size());
		MPI_Waitall((int)vReceiveRequests.size(), vReceiveRequests.data(), vReceiveStatuses.data());

		return vSendRequests;
	}

	MPI_Datatype cgm_mpi::create_datatype(std::vector<int> indexes) const
	{
		MPI_Datatype res;
		size_t receive_size = indexes.size();
		std::vector<int> receive_lens(receive_size, 1);
		MPI_Type_indexed((int)receive_size, receive_lens.data(), indexes.data(), MPI_INT, &res);
		return res;
	}

	cgm_mpi::side_index cgm_mpi::get_opposite_side(cgm_mpi::side_index side) const
	{
		switch(side)
		{
			case side_index::LEFT:
				return side_index::RIGHT;
			case side_index::RIGHT:
				return side_index::LEFT;
			case side_index::TOP:
				return side_index::BOTTOM;
			case side_index::BOTTOM:
				return side_index::TOP;
			default:
				throw std::runtime_error("Wrong side index: " + std::to_string(side));
		}
	}
}

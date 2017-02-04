#include <mpi.h>
#include <sstream>
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
			this->vBoundaryProcs[i].nProcessRank = boundInfo.nProcessRank;
			if (boundInfo.nProcessRank == fdm_grid_mpi::NOT_PROCESS)
			{
				this->vBoundaryProcs[i].receive_t = nullptr;
				this->vBoundaryProcs[i].send_t = nullptr;
				continue;
			}

//			int nRank;
//			MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
//
//			std::cout << nRank << ": " << i << ": ";
//			for (int receive : boundInfo.vBoundariesToReceive)
//			{
//				std::cout << receive << " ";
//			}
//			std::cout << std::endl;

			this->create_datatype(boundInfo.vBoundariesToReceive, &this->vBoundaryProcs[i].receive_t);
			this->create_datatype(boundInfo.vBoundariesToSend, &this->vBoundaryProcs[i].send_t);
		}

		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		std::ostringstream str;
		if (nRank == 0)
		{
			for (auto &&val : vBoundaryProcs[side_index::TOP].vBoundariesToReceive)
			{
				str << val << "\t";
			}
			std::cout << "0 indexes: " << str.str() << std::endl;
		}

		if (nRank == 5)
		{
			for (auto &&val : vBoundaryProcs[side_index::BOTTOM].vBoundariesToSend)
			{
				str << val << "\t";
			}
			std::cout << "5 indexes: " << str.str() << std::endl;
		}
	}

	size_t cgm_mpi::solve(const sparse_matrix &mA, vector &vSol, const vector &vRightPart, double rEps)
	{
		vSol = 0.0;

		vector vGradient = vRightPart - mA.vector_multiply(vSol);
		vector vResiduals = vGradient;
		double rBRR1 = vGradient * vGradient;

		size_t nIterationNum;
		std::vector<MPI_Request> vSendRequests(this->vBoundaryProcs.size(), MPI_REQUEST_NULL);
		std::vector<MPI_Status> vSendStatuses(vSendRequests.size());

		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		if (nRank == 0)
		{
			std::cout << "UPPER: " << this->vBoundaryProcs[side_index::TOP].nProcessRank << std::endl;
		}

		if (nRank == 5)
		{
			std::cout << "LOWER: " << this->vBoundaryProcs[side_index::BOTTOM].nProcessRank << std::endl;
		}

		double rGradNorm = vGradient.norm_inf();
		double rGradNormFull = rGradNorm;
		for (nIterationNum = 0; (nIterationNum < 5) && (rGradNormFull > rEps); ++nIterationNum)
		{
//			std::cout << nRank << ": Iteration " << nIterationNum << "; Norm: " << rGradNorm << "; Norm full: " << rGradNormFull << std::endl;
			vector vTemp = mA.vector_multiply(vResiduals);
			double rAlpha = rBRR1 / (vTemp * vResiduals);

//			std::cout << nRank << ": Iteration " << nIterationNum << " 1" << std::endl;
			MPI_Waitall((int)vSendRequests.size(), vSendRequests.data(), vSendStatuses.data());
			vSol.add_scale_vector(vResiduals, rAlpha);

//			std::cout << nRank << ": Iteration " << nIterationNum << " 2" << std::endl;

			vSendRequests = this->update_boundaries(vSol);

//			std::cout << nRank << ": Iteration " << nIterationNum << " 3" << std::endl;
			vGradient = vRightPart - mA.vector_multiply(vSol); // TODO: Optimize (try to return simple implementation)

			double rBRR2  = vGradient * vGradient;
			double rBetta = rBRR2 / rBRR1;
			rBRR1  = rBRR2;

			vTemp = vGradient;
			vTemp.add_scale_vector(vResiduals, rBetta);
			vResiduals = vTemp;

			rGradNorm = vGradient.norm_inf();
			MPI_Allreduce(&rGradNorm, &rGradNormFull, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		}

		std::cout << nRank << " Iteration " << nIterationNum << " END!! Norm: " << rGradNorm << "; Norm full: " << rGradNormFull << std::endl;
		return nIterationNum;
	}

	std::vector<MPI_Request> cgm_mpi::update_boundaries(vector &vSolution) const
	{
		std::vector<MPI_Request> vSendRequests(this->vBoundaryProcs.size(), MPI_REQUEST_NULL);
		std::vector<MPI_Request> vReceiveRequests(this->vBoundaryProcs.size(), MPI_REQUEST_NULL);

		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
//		std::ostringstream str;
//		for (auto const &proc : this->vBoundaryProcs)
//		{
//			str << proc.nProcessRank << ' ';
//		}
//		std::cout << "Ranks (" << nRank << "): " << str.str() << std::endl;
		if (nRank == 5)
		{
//			std::cout << vSolution.size() << std::endl;
			int nRows = 27, nColumns = 11;
			std::ostringstream str;
			str << nRank << " (-1) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[i] << "\t";
			}
			str << std::endl;

			str << nRank << " (-2) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[nColumns + i] << "\t";
			}
			str << std::endl;
			std::cout << str.str() << std::flush;
		}

		if (nRank == 0)
		{
//			std::cout << vSolution.size() << std::endl;
			std::ostringstream str;
			int nRows = 26, nColumns = 11;
			str << nRank << " (1) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[(nRows - 1) * nColumns + i] << "\t";
//				str << (nRows - 1) * nColumns + i << "\t";
			}
			str << std::endl;

			str << nRank << " (2) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[(nRows - 2) * nColumns + i] << "\t";
			}
			str << std::endl;
			std::cout << str.str() << std::flush;
		}

		for (int i = 0; i < this->vBoundaryProcs.size(); ++i)
		{
			auto const &proc = this->vBoundaryProcs[i];
			if (proc.nProcessRank == fdm_grid_mpi::NOT_PROCESS)
				continue;

//			std::cout << nRank << ": receive_t " << proc.receive_t  << std::endl;
			MPI_Irecv(&vSolution[0], 1, proc.receive_t, proc.nProcessRank, 0, MPI_COMM_WORLD, &vReceiveRequests[i]);
		}

//		std::cout << nRank << ": Irecv completed."  << std::endl;
		for (int i = 0; i < this->vBoundaryProcs.size(); ++i)
		{
			auto const &proc = this->vBoundaryProcs[i];
			if (proc.nProcessRank == fdm_grid_mpi::NOT_PROCESS)
				continue;

			MPI_Isend(&vSolution[0], 1, proc.send_t, proc.nProcessRank, 0, MPI_COMM_WORLD, &vSendRequests[i]);
		}

//		std::cout << nRank << ": Isend completed."  << std::endl;

		std::vector<MPI_Status> vReceiveStatuses(vReceiveRequests.size());
		MPI_Waitall((int)vReceiveRequests.size(), vReceiveRequests.data(), vReceiveStatuses.data());
		MPI_Waitall((int)vSendRequests.size(), vSendRequests.data(), vReceiveStatuses.data());

		if (nRank == 5)
		{
			std::cout << std::endl << std::endl;
//			std::cout << vSolution.size() << std::endl;
			int nRows = 27, nColumns = 11;
			std::ostringstream str;
			str << nRank << " (-1) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[i] << "\t";
			}
			str << std::endl;

			str << nRank << " (-2) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[nColumns + i] << "\t";
			}
			str << std::endl;
			std::cout << str.str() << std::flush;
			std::cout << std::endl << std::endl;
		}

		if (nRank == 0)
		{
//			std::cout << vSolution.size() << std::endl;
			std::ostringstream str;
			int nRows = 26, nColumns = 11;
			str << nRank << " (1) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[(nRows - 1) * nColumns + i] << "\t";
			}
			str << std::endl;

			str << nRank << " (2) ";
			for (int i = 0; i < nColumns; ++i)
			{
				str << vSolution[(nRows - 2) * nColumns + i] << "\t";
			}
			str << std::endl;
			std::cout << str.str() << std::flush;

			std::cout << std::endl << std::endl;
		}

		return vSendRequests;
	}

	void cgm_mpi::create_datatype(std::vector<int> indexes, MPI_Datatype *datatype) const
	{
		std::vector<int> element_lens(indexes.size(), 1);
		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		MPI_Type_indexed((int)indexes.size(), element_lens.data(), indexes.data(), MPI_INT, datatype);
		MPI_Type_commit(datatype);
//		std::cout << nRank << ": " << indexes.size() << " " << element_lens.size() << "; res: " << res << std::endl;
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

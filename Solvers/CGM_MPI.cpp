#include <mpi.h>
#include <sstream>
#include "CGM_MPI.h"

#include "SparseMatrix.h"

namespace Solvers
{
	cgm_mpi::cgm_mpi(const std::vector<fdm_grid_mpi::boundaries_list> &vBoundaryProcs,
		                 const std::vector<std::vector<size_t>> &vNeighbours, const std::vector<bool> &vBoundFlag)
		: m_vBoundaryProcs(vBoundaryProcs.size())
		, m_vNeighbours(vNeighbours)
		, m_vBoundFlag(vBoundFlag)
	{
		for (int i = 0; i < side_index::SIDES_NUMBER; ++i)
		{
			auto const &boundInfo = vBoundaryProcs[i];
			this->m_vBoundaryProcs[i].nProcessRank = boundInfo.nProcessRank;
			if (boundInfo.nProcessRank == fdm_grid_mpi::NOT_PROCESS)
			{
				this->m_vBoundaryProcs[i].receive_t = nullptr;
				this->m_vBoundaryProcs[i].send_t = nullptr;
				continue;
			}

			this->create_datatype(boundInfo.vBoundariesToReceive, &this->m_vBoundaryProcs[i].receive_t);
			this->create_datatype(boundInfo.vBoundariesToSend, &this->m_vBoundaryProcs[i].send_t);
		}
	}

	size_t cgm_mpi::solve(const sparse_matrix &mA, vector &vSol, const vector &vRightPart, double rEps, int nMaxIter)
	{
		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		vSol = 0.0;

		vector vGradient = vRightPart - mA.vector_multiply(vSol);
		vector vResiduals = vGradient;
		double rBRR1 = vGradient * vGradient;

		size_t nIterationNum = 0;
		double rGradNorm = vGradient.norm_inf();
		double rGradNormFull = rGradNorm;
		while ((nIterationNum < nMaxIter) && (rGradNormFull > rEps))
		{
			vector vTemp = mA.vector_multiply(vResiduals);
			double rAlpha = rBRR1 / (vTemp * vResiduals);

			vSol.add_scale_vector(vResiduals, rAlpha);

			vector vRightPartDelta;
			this->update_boundaries(vSol, vRightPartDelta);
			vGradient = vRightPart + vRightPartDelta - mA.vector_multiply(vSol);

			double rBRR2  = vGradient * vGradient;
			double rBetta = rBRR2 / rBRR1;
			rBRR1  = rBRR2;

			vTemp = vGradient;
			vTemp.add_scale_vector(vResiduals, rBetta);
			vResiduals = vTemp;

			rGradNorm = vGradient.norm_inf();

			++nIterationNum;
			MPI_Allreduce(&rGradNorm, &rGradNormFull, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		}

		std::cout << nRank << " Iteration " << nIterationNum << " END. Norm: " << rGradNorm << "; Norm full: " << rGradNormFull << std::endl;
		return nIterationNum;
	}

	void cgm_mpi::update_boundaries(vector &vSolution, vector &vRightPartDelta) const
	{
		std::vector<MPI_Request> vSendRequests(this->m_vBoundaryProcs.size(), MPI_REQUEST_NULL);
		std::vector<MPI_Request> vReceiveRequests(this->m_vBoundaryProcs.size(), MPI_REQUEST_NULL);

		std::vector<double> vBoundaryValues(vSolution.size());

		for (int i = 0; i < this->m_vBoundaryProcs.size(); ++i)
		{
			auto const &proc = this->m_vBoundaryProcs[i];
			if (proc.nProcessRank == fdm_grid_mpi::NOT_PROCESS)
				continue;

			MPI_Irecv(vBoundaryValues.data(), 1, proc.receive_t, proc.nProcessRank, 0, MPI_COMM_WORLD, &vReceiveRequests[i]);
		}

		for (int i = 0; i < this->m_vBoundaryProcs.size(); ++i)
		{
			auto const &proc = this->m_vBoundaryProcs[i];
			if (proc.nProcessRank == fdm_grid_mpi::NOT_PROCESS)
				continue;

			MPI_Isend(&vSolution[0], 1, proc.send_t, proc.nProcessRank, 0, MPI_COMM_WORLD, &vSendRequests[i]);
		}


		std::vector<MPI_Status> vReceiveStatuses(vReceiveRequests.size()), vSendStatuses(vReceiveRequests.size());
		MPI_Waitall((int)vSendRequests.size(), vSendRequests.data(), vSendStatuses.data());
		MPI_Waitall((int)vReceiveRequests.size(), vReceiveRequests.data(), vReceiveStatuses.data());

		this->update_right_part(vBoundaryValues, vRightPartDelta);
	}

	void cgm_mpi::create_datatype(std::vector<int> indexes, MPI_Datatype *datatype) const
	{
		std::vector<int> element_lens(indexes.size(), 1);
		int nRank;
		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		MPI_Type_indexed((int)indexes.size(), element_lens.data(), indexes.data(), MPI_DOUBLE, datatype);
		MPI_Type_commit(datatype);
	}

	template <typename T>
	std::ostream& cgm_mpi::print_vector(const T &vArr, int nIndDisplace, int nMaxInd, const std::string &sPrefix,
	                                    bool bOutIndex, const std::string &sDelim) const
	{
		if (nMaxInd < 0)
		{
			nMaxInd = vArr.size() - nIndDisplace;
		}
		std::ostringstream str;
		if (sPrefix != "")
		{
			str << sPrefix << sDelim;
		}

		if (bOutIndex)
		{
			for (int i = 0; i < nMaxInd; ++i)
			{
				str << nIndDisplace + i << sDelim;
			}
		}
		else
		{
			for (int i = 0; i < nMaxInd; ++i)
			{
				str << vArr[nIndDisplace + i] << sDelim;
			}
		}

		return std::cout << str.str() << std::endl;
	}

	void cgm_mpi::update_right_part(const std::vector<double> &vBoundaryValues, vector &vRightPartDelta) const
	{
		vRightPartDelta.resize(vBoundaryValues.size(), 0);

		for (int nInd = 0; nInd < vBoundaryValues.size(); nInd++)
		{
			double rDelta = vBoundaryValues[nInd];
			if (std::abs(rDelta) < 1e-7)
				continue;

			vRightPartDelta[nInd] = rDelta;

			for (size_t nNeighbourInd : this->m_vNeighbours.at(nInd))
			{
				if (nNeighbourInd == fdm_grid_mpi::EMPTY_NEIGHBOUR || this->m_vBoundFlag.at(nNeighbourInd))
					continue;

				vRightPartDelta[nNeighbourInd] += rDelta;
			}
		}
	}
}

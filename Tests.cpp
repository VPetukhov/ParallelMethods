#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tests

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <Solvers/CGM_MPI.h>
#include <Utils/FDMGrid.h>
#include <Utils/TaskFunctions.h>

#include "Solvers/Functions.h"
#include "Utils/SparseMatrix.h"

BOOST_AUTO_TEST_SUITE(Tests)

	BOOST_AUTO_TEST_CASE(test_vector)
	{
		double d1[5] = {1, 2, 3, 4, 5};
		double d2[5] = {0, 1, 2, 3, 4};
		vector v1(std::valarray<double>(d1, 5));
		vector v2(std::valarray<double>(d2, v1.size()));


		BOOST_CHECK((v1 - v1).apply(std::fabs).max() < 1e-7);
		BOOST_CHECK(std::abs((v1 - v2).apply(std::fabs).max() - 1) < 1e-7);
	}

	BOOST_AUTO_TEST_CASE(test_sparse_matrix)
	{
		const int rank = 5;
		double a[rank][rank] = {
				{0.466004,-0.206244,-0.0322768,0.262886,0.13187},
				{-0.206244,0.439495,0.027306,-0.0494091,-0.28005},
				{-0.0322768,0.027306,0.269317,-0.00742683,-0.0431571},
				{0.262886,-0.0494091,-0.00742683,0.439013,0.12912},
				{0.13187,-0.28005,-0.0431571,0.12912,0.492918}
		};

		sparse_matrix mA(rank, rank);
		for (size_t i = 0; i < rank; ++i)
		{
			for (size_t j = 0; j < rank; ++j)
			{
				mA.at(i, j) = a[i][j];
			}
		}

		for (size_t i = 0; i < rank; ++i)
		{
			for (size_t j = 0; j < rank; ++j)
			{
				BOOST_CHECK(mA.at(i, j) == a[i][j]);
			}
		}
	}

	BOOST_AUTO_TEST_CASE(test_pcgm_trivial)
	{
		const int equationRank = 4;
		const double eps = 1e-7;

		sparse_matrix mA(equationRank, 3);
		vector vRightPart(equationRank);
		for (size_t i = 0; i < equationRank; ++i)
		{
			mA.at(i, i) = 1;
			vRightPart[i] = 1;
		}
		vRightPart[1] = 2;

		vector vSol(equationRank);
		Solvers::CGM(mA, vSol, vRightPart, eps);

		for (size_t i = 0; i < equationRank; ++i)
		{
			BOOST_CHECK(std::abs(vSol[i] - vRightPart[i]) < eps);
		}
	}

	BOOST_AUTO_TEST_CASE(test_pcgm)
	{
		const int equationRank = 5;
		const double eps = 1e-5;

		double a[equationRank][equationRank] = {
				{0.466004,-0.206244,-0.0322768,0.262886,0.13187},
				{-0.206244,0.439495,0.027306,-0.0494091,-0.28005},
				{-0.0322768,0.027306,0.269317,-0.00742683,-0.0431571},
				{0.262886,-0.0494091,-0.00742683,0.439013,0.12912},
				{0.13187,-0.28005,-0.0431571,0.12912,0.492918}
		};
		double b[equationRank] = {0.565811, 0.505768, 0.816686, 0.584653, 0.0253342};
		double solution[equationRank] = {2.44934048, 2.98386069, 3.24620609, -0.16122111, 1.41784955};
		vector vSolution(std::valarray<double>(solution, equationRank));


		sparse_matrix mA(equationRank, equationRank);
		vector vRightPart(equationRank);
		for (size_t i = 0; i < equationRank; ++i)
		{
			for (size_t j = 0; j < equationRank; ++j)
			{
				mA.at(i, j) = a[i][j];
			}
			vRightPart[i] = b[i];
		}

		vector vSol(equationRank);
		Solvers::CGM(mA, vSol, vRightPart, eps / 1000);

		BOOST_CHECK((vSol - vSolution).apply(std::fabs).max() < eps);
		BOOST_CHECK((mA.vector_multiply(vSol) - vRightPart).norm_inf() < eps);
	}

BOOST_AUTO_TEST_SUITE_END()
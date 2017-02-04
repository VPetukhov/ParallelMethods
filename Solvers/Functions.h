#pragma once

#include "SparseMatrix.h"
namespace Solvers
{
	extern size_t CGM(const sparse_matrix &mA, vector &vSol, const vector &vRightPart, double rEps);
}
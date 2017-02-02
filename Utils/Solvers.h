#pragma once

#include "SparseMatrix.h"

extern size_t PCGM(const sparse_matrix &mA, vector &vSol, const vector &vRightPart, double rEps);
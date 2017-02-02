#pragma once

#include "SparseMatrix.h"

extern void PCGM(const sparse_matrix &mA, vector &vSol, const vector &vRightPart, double rEps);
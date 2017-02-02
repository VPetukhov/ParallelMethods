#pragma once

#include "FDMGrid.h"

namespace Task
{
	extern double exact_solution(const fdm_grid::coord &c);
	extern double bound_condition(const fdm_grid::coord &c);
	extern double right_part(const fdm_grid::coord &c);
}
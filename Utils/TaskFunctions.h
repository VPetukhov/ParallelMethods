#pragma once

#include "FDMGrid.h"

namespace Task
{
	extern double exact_solution(const coord &c);
	extern double bound_condition(const coord &c);
	extern double right_part(const coord &c);
}
#include "TaskFunctions.h"

#include "cmath"

namespace Task
{
	double exact_solution(const fdm_grid::coord &c)
	{
		return sin(c.rX + c.rY);
	}

	double bound_condition(const fdm_grid::coord &c)
	{
		return exact_solution(c);
	}

	double right_part(const fdm_grid::coord &c)
	{
		return 2.0 * sin(c.rX + c.rY);
	}

}
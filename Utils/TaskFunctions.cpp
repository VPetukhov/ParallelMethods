#include "math.h"

#include "TaskFunctions.h"

namespace Task
{
	double exact_solution(const coord &c)
	{
		return sin(c.rX + c.rY);
	}

	double bound_condition(const coord &c)
	{
		return exact_solution(c);
	}

	double right_part(const coord &c)
	{
		return 2.0 * sin(c.rX + c.rY);
	}

}
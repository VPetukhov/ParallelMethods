#include "Solvers.h"

void PCGM(const sparse_matrix& mA, vector& vSol, const vector& vRightPart, double rEps)
{
	size_t  nSize = vRightPart.size();

	vector vGradient(vRightPart), vResiduals(nSize), vTemp(nSize);

	vSol = 0.0;
	vTemp = mA.vector_multiply(vSol);
	vGradient.add_scale_vector(vTemp, -1.0);

	vResiduals = vGradient;
	double rBRR1 = vResiduals * vGradient;

	for (size_t iteration_num = 0; (iteration_num < nSize) && (vGradient.norm_inf() > rEps); ++iteration_num)
	{
		vTemp = mA.vector_multiply(vResiduals);
		double rAlpha = rBRR1 / (vTemp * vResiduals);

		vSol.add_scale_vector(vResiduals, rAlpha);
		vGradient.add_scale_vector(vTemp, - rAlpha);

		double rBRR2  = vGradient * vGradient;
		double rBetta = rBRR2 / rBRR1;
		rBRR1  = rBRR2;

		vTemp = vGradient;
		vTemp.add_scale_vector(vResiduals, rBetta);
		vResiduals = vTemp;
	}
}
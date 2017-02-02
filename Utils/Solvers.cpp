#include "Solvers.h"

size_t PCGM(const sparse_matrix& mA, vector& vSol, const vector& vRightPart, double rEps)
{
	size_t  nSize = vRightPart.size();

	vector vGradient(vRightPart), vResiduals(nSize), vTemp(nSize);

	vSol = 0.0;
	vTemp = mA.vector_multiply(vSol);
	vGradient.add_scale_vector(vTemp, -1.0);

	vResiduals = vGradient;
	double rBRR1 = vResiduals * vGradient;

	size_t nIterationNum;
	for (nIterationNum = 0; (nIterationNum < nSize) && (vGradient.norm_inf() > rEps); ++nIterationNum)
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
	return nIterationNum;
}
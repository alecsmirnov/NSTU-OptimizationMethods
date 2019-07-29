#include "goldenRatioMethod.h"

#include "minIntervalSearch.h"

#include <math.h>

GRMResult goldenRatioMethod(double (*func)(double, double), const double x[GR_APPROACH_SIZE], 
							const double S[GR_APPROACH_SIZE], double lambda0, double eps) {
	MISResult interval = minIntervalSearch(func, x, S, lambda0, eps);

	const double ratio_a = (3 - sqrt(5)) / 2;
	const double ratio_b = (sqrt(5) - 1) / 2;

	double lambda1 = interval.a + ratio_a * (interval.b - interval.a);
	double lambda2 = interval.a + ratio_b * (interval.b - interval.a);

	double fx1 = func(x[0] + lambda1 * S[0], x[1] + lambda1 * S[1]);
	double fx2 = func(x[0] + lambda2 * S[0], x[1] + lambda2 * S[1]);

	uint32_t calcs = interval.calcs_count + 2;
	while (eps < fabs(interval.b - interval.a)) {
		if (fx1 < fx2) {
			interval.b = lambda2;

			lambda2 = lambda1;
			fx2 = fx1;

			lambda1 = interval.a + ratio_a * (interval.b - interval.a);
			fx1 = func(x[0] + lambda1 * S[0], x[1] + lambda1 * S[1]);
		}
		else {
			interval.a = lambda1;

			lambda1 = lambda2;
			fx1 = fx2;

			lambda2 = interval.a + ratio_b * (interval.b - interval.a);
			fx2 = func(x[0] + lambda2 * S[0], x[1] + lambda2 * S[1]);
		}

		++calcs;
	}

	return (GRMResult){(interval.a + interval.b) / 2, calcs};
}

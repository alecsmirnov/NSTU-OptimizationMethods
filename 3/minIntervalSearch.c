#include "minIntervalSearch.h"

MISResult minIntervalSearch(double (*func)(double, double), const double x[MI_APPROACH_SIZE], const double S[MI_APPROACH_SIZE], double lambda0, double delta) {
	double fx0 = func(x[0] + lambda0 * S[0], x[1] + lambda0 * S[1]);

	if (fx0 < func(x[0] + (lambda0 + delta) * S[0], x[1] + (lambda0 + delta) * S[1]))
		delta = -delta;

	double lambda1 = lambda0 + delta;
	double fx1 = func(x[0] + lambda1 * S[0], x[1] + lambda1 * S[1]);

	uint32_t calcs = 3;
	do {
		delta *= 2;

		lambda0 = lambda1;
		lambda1 += delta;

		fx0 = fx1;
		fx1 = func(x[0] + lambda1 * S[0], x[1] + lambda1 * S[1]);

		++calcs;
	} while (fx1 < fx0);

	return lambda1 < lambda0 ? (MISResult){lambda1, lambda0, calcs} : (MISResult){lambda0, lambda1, calcs};
}
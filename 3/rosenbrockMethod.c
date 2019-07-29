#include "rosenbrockMethod.h"

#include "goldenRatioMethod.h"

#include <stdio.h>
#include <memory.h>

RMResult rosenbrockMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps) {
	double x[APPROACH_SIZE];
	memcpy(x, x0, sizeof(double) * APPROACH_SIZE);

	double S[APPROACH_SIZE][APPROACH_SIZE] = {{1, 0}, {0, 1}};

	uint32_t iters = 0;
	uint32_t calcs = 1;
	double prev_x[APPROACH_SIZE];
	double fx_prev;
	double fx;
	do {
		memcpy(prev_x, x, sizeof(double) * APPROACH_SIZE);

		double lambda[APPROACH_SIZE] = {0};
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i) {
			GRMResult gold_ratio = goldenRatioMethod(func, x, S[i], lambda[i], eps);
			lambda[i] = gold_ratio.func_min;
			calcs += gold_ratio.calcs_count;

			for (uint8_t j = 0; j != APPROACH_SIZE; ++j)
				x[j] += S[i][j] * lambda[i];
		}

		double A[APPROACH_SIZE][APPROACH_SIZE] = {0};
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			A[0][i] = S[0][i] * lambda[0] + S[1][i] * lambda[1];

		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			A[1][i] = lambda[1] * S[1][i];

		if (fabs(lambda[0]) < fabs(lambda[1]))
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			A[1][i] = lambda[0] * S[0][i];

		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			S[0][i] = A[0][i] / arrayNorm(A[0]);

		double K = A[1][0] * S[0][0] + A[1][1] * S[0][1];
		double B[APPROACH_SIZE];
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			B[i] = A[1][i] - K * S[0][i];

		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			S[1][i] = B[i] / arrayNorm(B);

		fx_prev = func(prev_x[0], prev_x[1]);
		fx = func(x[0], x[1]);
		calcs += 2;
		
		++iters;
	} while (eps < fabs(fx_prev - fx) &&
			 eps < fabs(prev_x[0] - x[0]) && eps < fabs(prev_x[1] - x[1]));

	return (RMResult){x[0], x[1], fx, iters, calcs};
}
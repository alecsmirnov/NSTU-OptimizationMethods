#include "directSearch.h"
#include "plotter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ROSENBROCK "Results/Rosenbrock/"
#define BROYDEN	   "Results/Broyden/"

#define arrayCopy(dest, src) dest[0] = src[0]; dest[1] = src[1]
#define arrayNorm(x) sqrt(x[0] * x[0] + x[1] * x[1])

typedef struct IntervalResult {
	double a;
	double b;

	uint32_t calcs_count;
} IntervalResult;

typedef struct GoldenRatioResult {
	double func_min;

	uint32_t calcs_count;
} GoldenRatioResult;

typedef struct GradResult {
	double x0[APPROACH_SIZE];

	uint32_t calcs_count;
} GradResult;

static char* concatStr(const char* s1, const char* s2) {
	char* result = (char*)malloc(sizeof(char) * (strlen(s1) + strlen(s2) + 1));

	strcpy(result, s1);
	strcat(result, s2);

	return result;
}

static void writeTableIter(FILE* fp, uint32_t iter, const double prev_x[APPROACH_SIZE], 
						   const double x[APPROACH_SIZE], double fx_prev, double fx, double lambda) {
	fprintf(fp, "%u\t", iter);
	fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", x[0], x[1]);
	fprintf(fp, "%."PRINT_ACCURACY"lf\t", fx);
	fprintf(fp, "%."PRINT_ACCURACY"lf\t", lambda);
	fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", x[0] - prev_x[0], x[1] - prev_x[1]);
	fprintf(fp, "%."PRINT_ACCURACY"lf\t", fx - fx_prev);
}

static void writeRosenbrockTableIter(FILE* fp, uint32_t iter, const double prev_x[APPROACH_SIZE], const double x[APPROACH_SIZE],
									 double fx_prev, double fx, const double S[APPROACH_SIZE][APPROACH_SIZE], double lambda) {
	writeTableIter(fp, iter, prev_x, x, fx_prev, fx, lambda);
	fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\n", S[0][0], S[0][1], S[1][0], S[1][1]);
}

static void writeBroydenTableIter(FILE* fp, uint32_t iter, const double prev_x[APPROACH_SIZE], const double x[APPROACH_SIZE], double fx_prev, 
								  double fx, const double A[APPROACH_SIZE][APPROACH_SIZE], const double grad[APPROACH_SIZE], double lambda) {
	writeTableIter(fp, iter, prev_x, x, fx_prev, fx, lambda);
	fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", grad[0], grad[1]);
	fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\n", A[0][0], A[0][1], A[1][0], A[1][1]);
}

static IntervalResult findIntervalMin(double (*func)(double, double), const double x[APPROACH_SIZE], const double S[APPROACH_SIZE], double lambda0) {
	double fx0 = func(x[0] + lambda0 * S[0], x[1] + lambda0 * S[1]);

	double delta = 1E-8;
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

	return lambda1 < lambda0 ? (IntervalResult){lambda1, lambda0, calcs} : (IntervalResult){lambda0, lambda1, calcs};
}

static GoldenRatioResult goldenRatioMethod(double (*func)(double, double), const double x[APPROACH_SIZE], const double S[APPROACH_SIZE], double lambda0, double eps) {
	IntervalResult interval = findIntervalMin(func, x, S, lambda0);

	const double ratio_a = (3 - sqrt(5)) / 2;
	const double ratio_b = (sqrt(5) - 3) / 2;

	double lambda1 = interval.a + ratio_a * (interval.b - interval.a);
	double lambda2 = interval.b + ratio_b * (interval.b - interval.a);

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

			lambda2 = interval.b + ratio_b * (interval.b - interval.a);
			fx2 = func(x[0] + lambda2 * S[0], x[1] + lambda2 * S[1]);
		}

		++calcs;
	}

	return (GoldenRatioResult){(interval.a + interval.b) / 2, calcs};
}

DMResult rosenbrockMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename) {
	char* filename_table = concatStr(ROSENBROCK, filename);
	FILE* fp = fopen(filename_table, "w");

	addStep(x0[0], x0[1]);

	double x[APPROACH_SIZE];
	arrayCopy(x, x0);

	double S[APPROACH_SIZE][APPROACH_SIZE] = {{1, 0}, {0, 1}};

	uint32_t iters = 0;
	uint32_t calcs = 1;
	double prev_x[APPROACH_SIZE];
	double fx_prev;
	double fx;
	do {
		arrayCopy(prev_x, x);

		double lambda[APPROACH_SIZE] = {0};
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i) {
			GoldenRatioResult gold_ratio = goldenRatioMethod(func, x, S[i], lambda[i], eps);
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
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			S[1][i] = A[1][i] - K * S[0][i];

		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			S[1][i] = S[1][i] / arrayNorm(S[1]);

		fx_prev = func(prev_x[0], prev_x[1]);
		fx = func(x[0], x[1]);
		calcs += 2;

		writeRosenbrockTableIter(fp, iters, prev_x, x, fx_prev, fx,  S, lambda[0]);
		addStep(x[0], x[1]);
		
		++iters;
	} while (eps < fabs(fx_prev - fx) && 
			 eps < fabs(prev_x[0] - x[0]) && eps < fabs(prev_x[1] - x[1]));

	fclose(fp);
	free(filename_table);
	
	makePicture(func, filename, ROSENBROCK);
	clearData();

	return (DMResult){x[0], x[1], fx, iters, calcs};
}

static GradResult grad(double(*func)(double, double), const double x1[APPROACH_SIZE], double eps) {
	const double fx1 = func(x1[0], x1[1]);

	double x0[APPROACH_SIZE];
	double x[APPROACH_SIZE];
	arrayCopy(x, x1);
	for (uint8_t i = 0; i != APPROACH_SIZE; ++i) {
		x[i] += eps;
		x0[i] = (func(x[0], x[1]) - fx1) / eps;
		x[i] -= eps;
	}

	return (GradResult){x0[0], x0[1], APPROACH_SIZE};
}

DMResult broydenMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename) {
	char* filename_table = concatStr(BROYDEN, filename);
	FILE* fp_table = fopen(filename_table, "w");

	addStep(x0[0], x0[1]);

	double x[APPROACH_SIZE];
	arrayCopy(x, x0);

	double eta[APPROACH_SIZE][APPROACH_SIZE] = {{1, 0}, {0, 1}};

	GradResult grad_result = grad(func, x0, eps);

	uint32_t iters = 0;
	uint32_t calcs = grad_result.calcs_count;
	double prev_x[APPROACH_SIZE];
	double fx_prev;
	double fx;
	do {
		double eta_grad[APPROACH_SIZE] = {0};
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			for (uint8_t j = 0; j != APPROACH_SIZE; ++j)
				eta_grad[i] += eta[i][j] * grad_result.x0[j];

		GoldenRatioResult lambda = goldenRatioMethod(func, x, eta_grad, 0, eps);
		calcs += lambda.calcs_count;

		double dx[APPROACH_SIZE];
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			dx[i] = lambda.func_min * eta_grad[i];

		arrayCopy(prev_x, x);
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			x[i] += dx[i];

		GradResult temp_grad = grad_result;
		grad_result = grad(func, x, eps);
		calcs += grad_result.calcs_count;

		double dgrad[APPROACH_SIZE];
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			dgrad[i] = grad_result.x0[i] - temp_grad.x0[i];

		double temp[APPROACH_SIZE];
		arrayCopy(temp, dx);
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			for (uint8_t j = 0; j != APPROACH_SIZE; ++j)
				temp[i] -= eta[i][j] * dgrad[j];
	
		double deta[APPROACH_SIZE][APPROACH_SIZE];
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			for (uint8_t j = 0; j != APPROACH_SIZE; ++j)
				deta[i][j] = temp[i] * temp[j] / (temp[0] * dgrad[0] + temp[1] * dgrad[1]);

		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			for (uint8_t j = 0; j != APPROACH_SIZE; ++j)
				eta[i][j] += deta[i][j];

		fx_prev = func(prev_x[0], prev_x[1]);
		fx = func(x[0], x[1]);
		calcs += 2;

		writeBroydenTableIter(fp_table, iters, prev_x, x, fx_prev, fx, eta, grad_result.x0, lambda.func_min);
		addStep(x[0], x[1]);

		++iters;
	} while (eps < fabs(arrayNorm(prev_x) - arrayNorm(x)) && 
			 eps < fabs(fx_prev - fx));

	fclose(fp_table);
	free(filename_table);

	makePicture(func, filename, BROYDEN);
	clearData();

	return (DMResult){x[0], x[1], func(x[0], x[1]), iters, calcs};
}
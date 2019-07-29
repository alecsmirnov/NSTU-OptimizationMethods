#include "directSearch.h"

#include "goldenRatioMethod.h"
#include "gradient.h"
#include "concatStr.h"
#include "plotter.h"

#include <stdio.h>

#define ROSENBROCK "Results/Rosenbrock/"
#define BROYDEN	   "Results/Broyden/"

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

DMResult rosenbrockMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename) {
	char* filename_table = concatStr(ROSENBROCK, filename);
	FILE* fp = fopen(filename_table, "w");

	plotterAddStep(x0[0], x0[1]);

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

		writeRosenbrockTableIter(fp, iters, prev_x, x, fx_prev, fx,  S, lambda[0]);
		plotterAddStep(x[0], x[1]);
		
		++iters;
	} while (eps < fabs(fx_prev - fx) && 
			 eps < fabs(prev_x[0] - x[0]) && eps < fabs(prev_x[1] - x[1]));

	fclose(fp);
	free(filename_table);
	
	plotterMakePicture(func, filename, ROSENBROCK);
	plotterClearData();

	return (DMResult){x[0], x[1], fx, iters, calcs};
}

DMResult broydenMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename) {
	char* filename_table = concatStr(BROYDEN, filename);
	FILE* fp_table = fopen(filename_table, "w");

	plotterAddStep(x0[0], x0[1]);

	double x[APPROACH_SIZE];
	memcpy(x, x0, sizeof(double) * APPROACH_SIZE);

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

		GRMResult lambda = goldenRatioMethod(func, x, eta_grad, 0, eps);
		calcs += lambda.calcs_count;

		double dx[APPROACH_SIZE];
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			dx[i] = lambda.func_min * eta_grad[i];

		memcpy(prev_x, x, sizeof(double) * APPROACH_SIZE);
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			x[i] += dx[i];

		GradResult temp_grad = grad_result;
		grad_result = grad(func, x, eps);
		calcs += grad_result.calcs_count;

		double dgrad[APPROACH_SIZE];
		for (uint8_t i = 0; i != APPROACH_SIZE; ++i)
			dgrad[i] = grad_result.x0[i] - temp_grad.x0[i];

		double temp[APPROACH_SIZE];
		memcpy(temp, dx, sizeof(double) * APPROACH_SIZE);
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
		plotterAddStep(x[0], x[1]);

		++iters;
	} while (eps < fabs(arrayNorm(prev_x) - arrayNorm(x)) && 
			 eps < fabs(fx_prev - fx));

	fclose(fp_table);
	free(filename_table);

	plotterMakePicture(func, filename, BROYDEN);
	plotterClearData();

	return (DMResult){x[0], x[1], func(x[0], x[1]), iters, calcs};
}
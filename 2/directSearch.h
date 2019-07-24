#ifndef DIRECTSEARCH_H
#define DIRECTSEARCH_H

#include <stdint.h>
#include <math.h>

#define APPROACH_SIZE 2
#define PRINT_ACCURACY "7"

typedef struct DMResult {
	double x[APPROACH_SIZE];
	double func_min;
	
	uint32_t iters_count;
	uint32_t calcs_count;
} DMResult;

inline double arrayNorm(const double x[APPROACH_SIZE]) {
	return sqrt(x[0] * x[0] + x[1] * x[1]);
}

DMResult rosenbrockMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename);
DMResult broydenMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename);

#endif
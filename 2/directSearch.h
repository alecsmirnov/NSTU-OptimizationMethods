#ifndef DIRECTSEARCH_H
#define DIRECTSEARCH_H

#include <stdint.h>

#define APPROACH_SIZE 2

#define PRECISION "7"

typedef struct DMResult {
	double x[APPROACH_SIZE];
	double func_min;
	
	uint32_t iters_count;
	uint32_t calcs_count;
} DMResult;

DMResult rosenbrockMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename);
DMResult broydenMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps, const char* filename);

#endif
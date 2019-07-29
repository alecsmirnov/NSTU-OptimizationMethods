#ifndef ROSENBROCKMETHOD_H
#define ROSENBROCKMETHOD_H

#include <stdint.h>
#include <math.h>

#include "plotter.h"

#define APPROACH_SIZE 2

typedef struct RMResult {
	double x[APPROACH_SIZE];
	double func_min;
	
	uint32_t iters_count;
	uint32_t calcs_count;
} RMResult;

inline double arrayNorm(const double x[APPROACH_SIZE]) {
	return sqrt(x[0] * x[0] + x[1] * x[1]);
}

RMResult rosenbrockMethod(double (*func)(double, double), const double x0[APPROACH_SIZE], double eps);

#endif
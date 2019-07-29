#ifndef GOLDENRATIOMETHOD_H
#define GOLDENRATIOMETHOD_H

#include <stdint.h>

#define GR_APPROACH_SIZE 2

typedef struct GRMResult {
	double func_min;

	uint32_t calcs_count;
} GRMResult;

GRMResult goldenRatioMethod(double (*func)(double, double), const double x[GR_APPROACH_SIZE], const double S[GR_APPROACH_SIZE], double lambda0, double eps);

#endif
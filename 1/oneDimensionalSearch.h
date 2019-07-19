#ifndef ONEDIMENSIONALSEARCH_H
#define ONEDIMENSIONALSEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct ODResult {
	double func_min;

	uint32_t iters_count;
	uint32_t calcs_count;
} ODResult;

typedef struct IntervalResult {
	double a;
	double b;
} IntervalResult;

ODResult dichotomyMethod(double (*func)(double), double a, double b, double eps, const char* filename);
ODResult goldenRatioMethod(double (*func)(double), double a, double b, double eps, const char* filename);
ODResult fibonacciMethod(double (*func)(double), double a, double b, double eps, const char* filename);
IntervalResult findIntervalMin(double (*func)(double), double x0, double delta, const char* filename);

#endif
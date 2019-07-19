#ifndef ONEDIMENSIONALSEARCH_H
#define ONEDIMENSIONALSEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct ODMResult {
	double func_min;

	uint32_t iters_count;
	uint32_t calcs_count;
} ODMResult;

typedef struct IntervalResult {
	double a;
	double b;
} IntervalResult;

ODMResult dichotomyMethod(double (*func)(double), double a, double b, double eps, const char* filename);
ODMResult goldenRatioMethod(double (*func)(double), double a, double b, double eps, const char* filename);
ODMResult fibonacciMethod(double (*func)(double), double a, double b, double eps, const char* filename);
IntervalResult findIntervalMin(double (*func)(double), double x0, double delta, const char* filename);

#endif
#ifndef ONEDIMENSIONALSEARCHMETHODS_H
#define ONEDIMENSIONALSEARCHMETHODS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct MethodResult {
	double func_min;

	uint32_t iters_count;
	uint32_t calcs_count;
} MethodResult;

typedef struct IntervalResult {
	double a;
	double b;
} IntervalResult;

MethodResult dichotomyMethod(double (*func)(double), double a, double b, double eps, const char* filename);
MethodResult goldenRatioMethod(double (*func)(double), double a, double b, double eps, const char* filename);
MethodResult fibonacciMethod(double (*func)(double), double a, double b, double eps, const char* filename);

IntervalResult findIntervalMin(double (*func)(double), double x0, double delta, const char* filename);

#endif
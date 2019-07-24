#ifndef ONEDIMENSIONALSEARCH_H
#define ONEDIMENSIONALSEARCH_H

#include <stdint.h>
#include <math.h>

typedef struct ODMResult {
	double func_min;

	uint32_t iters_count;
	uint32_t calcs_count;
} ODMResult;

typedef struct IntervalResult {
	double a;
	double b;
} IntervalResult;

inline double binet(uint32_t n) {
	const double phi = (1 + sqrt(5)) * 0.5;

	return round((pow(phi, n) - pow(1 - phi, n)) / sqrt(5));
}

ODMResult dichotomyMethod(double (*func)(double), double a, double b, double eps, const char* filename);
ODMResult goldenRatioMethod(double (*func)(double), double a, double b, double eps, const char* filename);
ODMResult fibonacciMethod(double (*func)(double), double a, double b, double eps, const char* filename);
IntervalResult findIntervalMin(double (*func)(double), double x0, double delta, const char* filename);

#endif
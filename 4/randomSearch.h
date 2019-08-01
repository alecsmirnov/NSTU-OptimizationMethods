#ifndef RANDOMSEARCH_H
#define RANDOMSEARCH_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct Point {
	double x;
	double y;
} Point;

typedef struct RSResult {
	double func_min;
	Point point_min;

	uint32_t iters_count;
	uint32_t calcs_count;
} RSResult;

inline double randInRange(double min, double max) {
	return min + rand() / (RAND_MAX / (max - min));
}

inline bool isInRange(double x0, double x, double x1) {
	return x0 < x && x < x1;
}

RSResult simpleRandomSearch(double (*func)(double, double), Point a, Point b, double P, double eps);

RSResult algorithm1(double (*func)(double, double), Point a, Point b, uint32_t m, double eps);
RSResult algorithm2(double (*func)(double, double), Point a, Point b, uint32_t m, double eps);
RSResult algorithm3(double (*func)(double, double), Point a, Point b, uint32_t m, double eps);

#endif
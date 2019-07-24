#ifndef MININTERVALSEARCH_H
#define MININTERVALSEARCH_H

#include <stdint.h>

#define MI_APPROACH_SIZE 2

typedef struct MISResult {
	double a;
	double b;

	uint32_t calcs_count;
} MISResult;

MISResult minIntervalSearch(double (*func)(double, double), const double x[MI_APPROACH_SIZE], const double S[MI_APPROACH_SIZE], double lambda0);

#endif

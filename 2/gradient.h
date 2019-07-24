#ifndef GRADIENT_H
#define GRADIENT_H

#include <stdint.h>
#include <memory.h>

#define VECTOR_SIZE 2

typedef struct GradResult {
	double x0[VECTOR_SIZE];

	uint32_t calcs_count;
} GradResult;

inline GradResult grad(double (*func)(double, double), const double x1[VECTOR_SIZE], double eps) {
	const double fx1 = func(x1[0], x1[1]);

	double x0[VECTOR_SIZE];
	double x[VECTOR_SIZE];
	memcpy(x, x1, sizeof(double) * VECTOR_SIZE);
	for (uint8_t i = 0; i != VECTOR_SIZE; ++i) {
		x[i] += eps;
		x0[i] = (func(x[0], x[1]) - fx1) / eps;
		x[i] -= eps;
	}

	return (GradResult){x0[0], x0[1], VECTOR_SIZE + 1};
}

#endif

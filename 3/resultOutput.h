#ifndef RESULTOUTPUT_H
#define RESULTOUTPUT_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#define PRINT_ACCURACY "7"
#define NONE 1

enum PrintColumnCode {
	PCC_ITER   = 1,
	PCC_CALC   = 2,
	PCC_X0	   = 4,
	PCC_X	   = 8,
	PCC_F	   = 16,
	PCC_R	   = 32,
	PCC_R_STEP = 64,
	PCC_ALPHA  = 128,
	PCC_EPS	   = 256,
	PCC_ALL    = 511
};

inline bool bitCheck(int val, uint8_t bit_pos) {
	return (val & (1 << bit_pos));
}

inline void writeTableHeader(FILE* fp, int val) {
	if (bitCheck(val, 0)) fprintf(fp, "i\t");
	if (bitCheck(val, 1)) fprintf(fp, "calcs\t");
	if (bitCheck(val, 2)) fprintf(fp, "x0\t\t\t");
	if (bitCheck(val, 3)) fprintf(fp, "x\t\t\t");
	if (bitCheck(val, 4)) fprintf(fp, "f(x)\t\t");
	if (bitCheck(val, 5)) fprintf(fp, "r\t\t");
	if (bitCheck(val, 6)) fprintf(fp, "r_step\t\t");
	if (bitCheck(val, 7)) fprintf(fp, "alpha\t");
	if (bitCheck(val, 8)) fprintf(fp, "eps");

	fprintf(fp, "\n");
}

inline void writeTableIter(FILE* fp, uint32_t iters, uint32_t calcs, const double x0[],
						   const double x[], double func_min, double r, double r_step, double eps, uint8_t alpha, int val) {
	if (bitCheck(val, 0)) fprintf(fp, "%u\t", iters);
	if (bitCheck(val, 1)) fprintf(fp, "%u\t", calcs);
	if (bitCheck(val, 2)) fprintf(fp, "%." PRINT_ACCURACY "lf %." PRINT_ACCURACY "lf\t", x0[0], x0[1]);
	if (bitCheck(val, 3)) fprintf(fp, "%." PRINT_ACCURACY "lf %." PRINT_ACCURACY "lf\t", x[0], x[1]);
	if (bitCheck(val, 4)) fprintf(fp, "%." PRINT_ACCURACY "e\t", func_min);
	if (bitCheck(val, 5)) fprintf(fp, "%." PRINT_ACCURACY "lf\t", r);
	if (bitCheck(val, 6)) fprintf(fp, "%." PRINT_ACCURACY "lf\t", r_step);
	if (bitCheck(val, 7)) fprintf(fp, "%u\t", alpha);
	if (bitCheck(val, 8)) fprintf(fp, "%.e", eps);

	fprintf(fp, "\n");
}

#endif
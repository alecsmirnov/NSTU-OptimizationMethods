#include <stdio.h>
#include <stdbool.h>
#include <memory.h>

#include "rosenbrockMethod.h"

#include "minIntervalSearch.h"
#include "goldenRatioMethod.h"

#define PRINT_ACCURACY "7"
#define SOLUTION_SEARCH_ITER_MAX 100

#define NONE 1

typedef double (*func_ptr)(double, double);

static struct Minimization {
	func_ptr G;
	double r;
	uint8_t alpha;
} minimization = {NULL, 0, 0};

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
	PCC_ALL	   = 511
};

static double f(double x, double y) {
	return pow(y - x, 2) + 10 * pow(x + 5, 2);
}

static double g(double x, double y) {
	return -(x + y);
}

// Penalty restriction functions
static double G1(double x, double y) {
	return 0 <= x + y ? 0 : -(x + y);					// a) x + y >= 0
}

static double G2(double x, double y) {			
	return pow(fabs(x + y - 1), minimization.alpha);	// b) x = 1 - y
}

// Barrier restriction functions
static double G3(double x, double y) {
	return g(x, y) <= 0 ? -1 / g(x, y) : INFINITY;
}

static double G4(double x, double y) {
	return g(x, y) <= 0 ? -log(-g(x, y)) : INFINITY;
}

static double G5(double x, double y) {
	return pow(0.5 * (g(x, y) + fabs(g(x, y))), minimization.alpha);
}

// Minimized function
static double Q(double x, double y) {
	return f(x, y) + minimization.r * minimization.G(x, y);
}

static void addRMResult(RMResult* dest, RMResult src, uint32_t iter) {
	memcpy(dest->x, src.x, sizeof(double) * APPROACH_SIZE);
	dest->calcs_count += src.calcs_count;
	dest->iters_count = iter;
	dest->func_min = src.func_min;
}

static RMResult restrictionRosenbrock(func_ptr G, double r, double r_step, uint8_t alpha, const double x0[APPROACH_SIZE], double eps) {
	RMResult result = {0, 0, 0, 0, 0};

	minimization.G = G;
	minimization.r = r;
	minimization.alpha = alpha;

	bool solution_found = false;
	for (uint32_t i = 0; i != SOLUTION_SEARCH_ITER_MAX && !solution_found; ++i) {
		addRMResult(&result, rosenbrockMethod(Q, x0, eps), i);
			
		if (minimization.r * minimization.G(result.x[0], result.x[1]) < eps)
			solution_found = true;
		else
			minimization.r *= r_step;
	}

	return result;
}

static bool bitCheck(int val, uint8_t bit_pos) {
	return ((val) & (1 << (bit_pos)));
}

static void writeTableHeader(FILE* fp, int val) {
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

static void writeTableIter(FILE* fp, uint64_t iters, uint64_t calcs, const double x0[APPROACH_SIZE], 
						   const double x[APPROACH_SIZE], double func_min, double r, double r_step, double eps, uint8_t alpha, int val) {	
	if (bitCheck(val, 0)) fprintf(fp, "%llu\t", iters);
	if (bitCheck(val, 1)) fprintf(fp, "%llu\t", calcs);
	if (bitCheck(val, 2)) fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", x0[0], x0[1]);
	if (bitCheck(val, 3)) fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", x[0], x[1]);
	if (bitCheck(val, 4)) fprintf(fp, "%."PRINT_ACCURACY"e\t", func_min);
	if (bitCheck(val, 5)) fprintf(fp, "%."PRINT_ACCURACY"lf\t", r);
	if (bitCheck(val, 6)) fprintf(fp, "%."PRINT_ACCURACY"lf\t", r_step);
	if (bitCheck(val, 7)) fprintf(fp, "%u\t", alpha);
	if (bitCheck(val, 8)) fprintf(fp, "%.e", eps);

	fprintf(fp, "\n");
}

static void researchEps(FILE* fp, func_ptr G, const double x0[APPROACH_SIZE], 
						double r, double r_step, double eps0, double eps1, double eps_step, int val) {
	writeTableHeader(fp, val);
	for (double eps = eps0; eps > eps1; eps *= eps_step) {
		RMResult result = restrictionRosenbrock(G, r, r_step, NONE, x0, eps);
		writeTableIter(fp, result.iters_count, result.calcs_count, x0, result.x, result.func_min, r, r_step, eps, NONE, val);
	}
}

static void researchX0(FILE* fp, func_ptr G, const double x0[][APPROACH_SIZE], uint8_t x0_size, double r, double r_step, int val) {
	const double eps = 1E-7;

	writeTableHeader(fp, val);
	for (uint8_t i = 0; i != x0_size; ++i) {
		RMResult result = restrictionRosenbrock(G, r, r_step, NONE, x0[i], eps);
		writeTableIter(fp, result.iters_count, result.calcs_count, x0[i], result.x, result.func_min, r, r_step, eps, NONE, val);
	}
}

static void researchRStep(FILE* fp, func_ptr G, const double x0[APPROACH_SIZE], double r, double r_step[], uint8_t r_step_size, int val) {
	const double eps = 1E-7;

	writeTableHeader(fp, val);
	for (uint8_t i = 0; i != r_step_size; ++i) {
		RMResult result = restrictionRosenbrock(G, r, r_step[i], NONE, x0, eps);
		writeTableIter(fp, result.iters_count, result.calcs_count, x0, result.x, result.func_min, r, r_step[i], eps, NONE, val);
	}
}

static void researchR(FILE* fp, func_ptr G, const double x0[APPROACH_SIZE], double r[], uint8_t r_size, double r_step, int val) {
	const double eps = 1E-7;

	writeTableHeader(fp, val);
	for (uint8_t i = 0; i != r_size; ++i) {
		RMResult result = restrictionRosenbrock(G, r[i], r_step, NONE, x0, eps);
		writeTableIter(fp, result.iters_count, result.calcs_count, x0, result.x, result.func_min, r[i], 5, eps, NONE, val);
	}
}

static void researchBarrierFuncAlpha(FILE* fp, const double x0[APPROACH_SIZE], double r, double r_step, uint8_t n_alpha, int val) {
	const double eps = 1E-7;

	writeTableHeader(fp, val);
	for (uint8_t i = 0; i != n_alpha; ++i) {
		RMResult result = restrictionRosenbrock(G5, 100000, 0.5, 2 * (i + 1), x0, eps);
		writeTableIter(stdout, result.iters_count, result.calcs_count, x0, result.x, result.func_min, r, r_step, eps, 2 * (i + 1), val);
	}
}

int main(int argc, char* argv[]) {
	double x0[APPROACH_SIZE] = {0, 1};
	double eps = 1E-7;

	printf("Unconditional task:\n");
	RMResult result = rosenbrockMethod(f, x0, eps);
	writeTableHeader(stdout, PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	writeTableIter(stdout, result.iters_count, result.calcs_count, x0, result.x, result.func_min, 0, 0, eps, NONE, PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("--------------------------------------------------------------------------------------------------------------");

	double r = 3;
	double r_step = 10;
	double eps0 = 1.0E-3;
	double eps1 = 1.0E-7;
	double eps_step = 1.0E-1;
	researchEps(stdout, G1, x0, r, r_step, eps0, eps1, eps_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 3;
	r_step = 10;
	enum {X0_PACK_SIZE = 9};
	double x0_pack[X0_PACK_SIZE][APPROACH_SIZE] = {{0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}, {2, 3}};
	researchX0(stdout, G1, x0_pack, X0_PACK_SIZE, r, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 5;
	enum {R_STEP_SIZE = 9};
	double r_step_pack[R_STEP_SIZE] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	researchRStep(stdout, G1, x0, r, r_step_pack, R_STEP_SIZE, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R_STEP | PCC_EPS);
	puts("");
	
	r_step = 6;
	enum {R_PACK_SIZE = 9};
	double r_pack[R_PACK_SIZE] = {2, 4, 6, 8, 10, 12, 14, 16, 18};
	researchR(stdout, G1, x0, r_pack, R_PACK_SIZE, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R | PCC_EPS);
	puts("--------------------------------------------------------------------------------------------------------------");

	r = 2;
	r_step = 4;
	eps0 = 1.0E-3;
	eps1 = 1.0E-7;
	eps_step = 1.0E-1;
	researchEps(stdout, G2, x0, r, r_step, eps0, eps1, eps_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 3;
	r_step = 10;
	//double x0_pack[X0_PACK_SIZE][APPROACH_SIZE] = {{0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}, {2, 3}};
	researchX0(stdout, G2, x0_pack, X0_PACK_SIZE, r, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 5;
	//double r_step_pack[R_STEP_SIZE] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	researchRStep(stdout, G2, x0, r, r_step_pack, R_STEP_SIZE, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R_STEP | PCC_EPS);
	puts("");

	r_step = 6;
	//double r_pack[R_PACK_SIZE] = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	researchR(stdout, G2, x0, r_pack, R_PACK_SIZE, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R | PCC_EPS);
	puts("--------------------------------------------------------------------------------------------------------------");

	r = 100;
	r_step = 0.5;
	eps0 = 1.0E-3;
	eps1 = 1.0E-7;
	eps_step = 1.0E-1;
	researchEps(stdout, G3, x0, r, r_step, eps0, eps1, eps_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 100;
	r_step = 0.5;
	//double x0_pack[X0_PACK_SIZE][APPROACH_SIZE] = {{0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}, {2, 3}};
	researchX0(stdout, G3, x0_pack, X0_PACK_SIZE, r, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 100;
	double r_step_pack2[R_STEP_SIZE] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
	researchRStep(stdout, G3, x0, r, r_step_pack2, R_STEP_SIZE, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R_STEP | PCC_EPS);
	puts("");

	r_step = 0.005;
	//double r_pack[R_PACK_SIZE] = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	researchR(stdout, G3, x0, r_pack, R_PACK_SIZE, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R | PCC_EPS);
	puts("--------------------------------------------------------------------------------------------------------------");

	r = 2;
	r_step = 4;
	eps0 = 1.0E-3;
	eps1 = 1.0E-7;
	eps_step = 1.0E-1;
	researchEps(stdout, G4, x0, r, r_step, eps0, eps1, eps_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 3;
	r_step = 10;
	//double x0_pack[X0_PACK_SIZE][APPROACH_SIZE] = {{0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}, {2, 3}};
	researchX0(stdout, G4, x0_pack, X0_PACK_SIZE, r, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_EPS);
	puts("");

	r = 3;
	//double r_step_pack2[R_STEP_SIZE] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
	researchRStep(stdout, G4, x0, r, r_step_pack2, R_STEP_SIZE, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R_STEP | PCC_EPS);
	puts("");

	r_step = 6;
	//double r_pack[R_PACK_SIZE] = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
	researchR(stdout, G4, x0, r_pack, R_PACK_SIZE, r_step, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_R | PCC_EPS);
	puts("--------------------------------------------------------------------------------------------------------------");

	r = 10000;
	r_step = 0.5;
	uint8_t n_alpha = 7;
	researchBarrierFuncAlpha(stdout, x0, r, r_step, n_alpha, PCC_ITER | PCC_CALC | PCC_X0 | PCC_X | PCC_F | PCC_ALPHA | PCC_EPS);
	puts("");

	return 0;
}
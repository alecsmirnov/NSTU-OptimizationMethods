#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include "randomSearch.h"

#define PRINT_ACCURACY "7"

typedef RSResult (*alg_ptr)(double (*func)(double, double), Point a, Point b, uint32_t m, double eps);

static double f(double x, double y) {
	enum {CONST_COUNT = 6};

	double C[CONST_COUNT] = { 6,  2,  4,  2,  8,  8};
	double a[CONST_COUNT] = {-3,  4, -8, -6,  3, -6};
	double b[CONST_COUNT] = { 9, -7,  3, -9, -2, -8};

	double result = 0;
	for (uint8_t i = 0; i != CONST_COUNT; ++i)
		result += C[i] / (1 + (x - a[i]) * (x - a[i]) + (y - b[i]) * (y - b[i]));

	return result;
}

static void simpleRandomSearchTable(Point a, Point b, const char* filename) {
	enum {EPS_PACK_SIZE = 2};
	double eps_pack[EPS_PACK_SIZE] = {1E-1, 1E-2};

	enum {P_PACK_SIZE = 9};
	double P_pack[P_PACK_SIZE] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

	FILE* fp = fopen(filename, "w");

	fprintf(fp, "eps\tP\tx\tfx\n");
	for (uint8_t i = 0; i != EPS_PACK_SIZE; ++i)
		for (uint8_t j = 0; j != P_PACK_SIZE; ++j) {
			RSResult result = simpleRandomSearch(f, a, b, P_pack[j], eps_pack[i]);

			fprintf(fp, "%.e\t", eps_pack[i]);
			fprintf(fp, "%."PRINT_ACCURACY"lf\t", P_pack[j]);
			fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", result.point_min.x, result.point_min.y);
			fprintf(fp, "%."PRINT_ACCURACY"lf\n", result.func_min);
		}

	fclose(fp);
}

static void algorithmTable(alg_ptr alg, Point a, Point b, double eps, const char* filename) {
	enum {M_PACK_SIZE = 5};
	uint32_t m_pack[M_PACK_SIZE] = {5, 10, 100, 500, 1000};

	FILE* fp = fopen(filename, "w");

	fprintf(fp, "m\tf_calcs\tx\tfx\n");
	for (uint8_t i = 0; i != M_PACK_SIZE; ++i) {
		RSResult result = alg(f, a, b, m_pack[i], eps);

		fprintf(fp, "%u\t", m_pack[i]);
		fprintf(fp, "%u\t", result.calcs_count);
		fprintf(fp, "%."PRINT_ACCURACY"e %."PRINT_ACCURACY"e\t", result.point_min.x, result.point_min.y);
		fprintf(fp, "%."PRINT_ACCURACY"lf\n", result.func_min);
	}

	fclose(fp);
}

int main(int argc, char* argv[]) {
	srand((uint32_t)time(NULL));

	Point a = {-10, -10};
	Point b = { 10,  10};
	double eps = 1E-3;

	simpleRandomSearchTable(a, b, "Results/simpleRandomSearch.txt");

	algorithmTable(algorithm1, a, b, eps, "Results/algorithm1.txt");
	algorithmTable(algorithm2, a, b, eps, "Results/algorithm2.txt");
	algorithmTable(algorithm3, a, b, eps, "Results/algorithm3.txt");

	return 0;
}
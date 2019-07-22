#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "directSearch.h"
#include "plotter.h"

static double f1(double x, double y) {
	double A1 = 1; 
	double A2 = 2;
	double a1 = 3;
	double a2 = 2;
	double b1 = 1;
	double b2 = 2;
	double c1 = 1;
	double c2 = 2;
	double d1 = 3;
	double d2 = 1;

	return  -(A1 / (1 + pow((x - a1) / b1, 2) + pow((y - c1) / d1, 2)) + 
			  A2 / (1 + pow((x - a2) / b2, 2) + pow((y - c2) / d2, 2)));
}

static double f2(double x, double y) {
	return 100 * pow(y - x, 2) + pow(1 - x, 2);
}

static double f3(double x, double y) {
	return 100 * pow(y - x * x, 2) + pow(1 - x, 2);
}

static writeGeneralTableIter(FILE* fp, double eps, uint8_t iters, uint8_t calcs, 
							 double x0[APPROACH_SIZE], double x[APPROACH_SIZE], double func_min) {
	fprintf(fp, "%."PRINT_ACCURACY"lf\t", eps);
	fprintf(fp, "%u\t", iters);
	fprintf(fp, "%u\t", calcs);
	fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", x0[0], x0[1]);
	fprintf(fp, "%."PRINT_ACCURACY"lf %."PRINT_ACCURACY"lf\t", x[0], x[1]);
	fprintf(fp, "%."PRINT_ACCURACY"lf\n", func_min);
}

static void makeTable(double (*func)(double, double), double x0[APPROACH_SIZE], const char* filename) {
	char gen_table[FILENAME_MAX];

	sprintf(gen_table, "Results/Rosenbrock/General_%s", filename);
	FILE* fp_rosenbrock = fopen(gen_table, "w");

	sprintf(gen_table, "Results/Broyden/General_%s", filename);
	FILE* fp_broyden = fopen(gen_table, "w");

	for (double eps = 1E-3; eps >= 1E-7; eps /= 10) {
		char table[FILENAME_MAX];
		sprintf(table, "Rosenbrock_%.0E_%s", eps, filename);

		DMResult result = rosenbrockMethod(func, x0, eps, table);
		writeGeneralTableIter(fp_rosenbrock, eps, result.iters_count, result.calcs_count, x0, result.x, result.func_min);

		sprintf(table, "Broyden_%.0E_%s", eps, filename);

		result = broydenMethod(func, x0, eps, table);
		writeGeneralTableIter(fp_broyden, eps, result.iters_count, result.calcs_count, x0, result.x, result.func_min);
	}

	fclose(fp_rosenbrock);
	fclose(fp_broyden);
}

int main(int argc, char* argv[]) {
	double x0[APPROACH_SIZE] = {1, 1};
	double eps = 1E-7;

	makeTable(f1, x0, "f1(1,1).txt");
	makeTable(f2, x0, "f2(1,1).txt");
	makeTable(f3, x0, "f3(1,1).txt");

	//rosenbrockMethod(f1, x0, eps, "Rosenbrock_fx1.txt");
	//rosenbrockMethod(f2, x0, eps, "Rosenbrock_fx2.txt");
	//rosenbrockMethod(f3, x0, eps, "Rosenbrock_fx3.txt");
	
	//broydenMethod(f1, x0, eps, "Broyden_fx1.txt");
	//broydenMethod(f2, x0, eps, "Broyden_fx2.txt");
	//broydenMethod(f3, x0, eps, "Broyden_fx3.txt");

	return 0;
}
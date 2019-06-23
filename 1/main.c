#include <stdio.h>
#include <math.h>

#include "oneDimensionalSearch.h"

static double func(double x) {
	return (x - 2) * (x - 2);
}

static void createLogTable(double a, double b, double x0, double eps0, double eps1, double eps_step) {
	FILE* fp = fopen("Results/logTable.txt", "w");
	fprintf(fp, "ln(eps)\tdich\tgoldRat\tfib\n");

	for (double eps = eps0; eps > eps1; eps *= eps_step) {
		MethodResult dich, gold, fib;
		char filename[FILENAME_MAX];

		sprintf(filename, "Results/dichotomy_%e.txt", eps);
		dich = dichotomyMethod(func, a, b, eps, filename);

		sprintf(filename, "Results/goldenRatio_%e.txt", eps);
		gold = goldenRatioMethod(func, a, b, eps, filename);

		sprintf(filename, "Results/fibonacci_%e.txt", eps);
		fib = fibonacciMethod(func, a, b, eps, filename);

		fprintf(fp, "%.14lf\t%u\t%u\t%u\n", log(eps), dich.calcs_count, gold.calcs_count, fib.calcs_count);
	}

	findIntervalMin(func, x0, eps0, "Results/findInterval.txt");

	fclose(fp);
}

int main(int argc, char* argv[]) {
	double a = -2;
	double b = 20;

	double x0 = -20;

	double eps0 = 1.0E-1;
	double eps1 = 1.0E-7;
	double eps_step = 1.0E-1;

	createLogTable(a, b, x0, eps0, eps1, eps_step);

	return 0;
}
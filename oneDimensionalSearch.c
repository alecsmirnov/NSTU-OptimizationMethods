#include "oneDimensionalSearch.h"

#include <math.h>

#define dichotomyCalculationNum(iter) ((iter) * 2)
#define goldenRatioCalculationNum(iter) ((iter) + 1)
#define fibonacciCalculationNum(iter) ((iter) + 1)

static FILE* openFile(const char* filename) {
	FILE* fp = fopen(filename, "w");

	if (fp == NULL) {
		fprintf(stderr, "Error: can't open file %s!\n", filename);
		exit(EXIT_FAILURE);
	}

	return fp;
}

static void closeFile(FILE* fp) {
	int close_status = fclose(fp);

	if (close_status == EOF) {
		fprintf(stderr, "Error: can't close file!\n");
		exit(EXIT_FAILURE);
	}
}

static double binet(uint32_t n) {
	const double phi = (1 + sqrt(5)) * 0.5;

	return round((pow(phi, n) - pow(1 - phi, n)) / sqrt(5));
}

MethodResult dichotomyMethod(double (*func)(double), double a, double b, double eps, const char* filename) {
	FILE* fp = openFile(filename);
	fprintf(fp, "i\ta\tb\t(b - a)\t(prev_b - prev_a)/(b - a)\tx1\tx2\tf(x1)\tf(x2)\n");

	const double delta = eps / 2;

	double prev_a = a;
	double prev_b = b;

	uint32_t iter = 0;
	while (eps < fabs(b - a)) {
		prev_a = a;
		prev_b = b;

		double x1 = (a + b - delta) / 2;
		double x2 = (a + b + delta) / 2;

		double fx1 = func(x1);
		double fx2 = func(x2);

		fx1 < fx2 ? (b = x2) : (a = x1);

		fprintf(fp, "%u\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\n",
				iter, a, b, b - a, (prev_b - prev_a) / (b - a), x1, x2, fx1, fx2);

		++iter;
	}

	fprintf(fp, "\neps:\t%.14lf\nxmin:\t%.14lf\ni:\t\t%u\nn:\t\t%u\n",
			eps, (a + b) / 2, iter, dichotomyCalculationNum(iter));

	closeFile(fp);

	return (MethodResult){(a + b) / 2, iter, dichotomyCalculationNum(iter)};
}

MethodResult goldenRatioMethod(double (*func)(double), double a, double b, double eps, const char* filename) {
	FILE* fp = openFile(filename);
	fprintf(fp, "i\ta\tb\t(b - a)\t(prev_b - prev_a)/(b - a)\tx1\tx2\tf(x1)\tf(x2)\n");

	const double ratio_a = (3 - sqrt(5)) / 2;
	const double ratio_b = (sqrt(5) - 3) / 2;

	double x1 = a + ratio_a * (b - a);
	double x2 = b + ratio_b * (b - a);

	double fx1 = func(x1);
	double fx2 = func(x2);

	uint32_t iter = 0;
	while (eps < fabs(b - a)) {
		double prev_a = a;
		double prev_b = b;

		if (fx1 < fx2) {
			b = x2;

			x2 = x1;
			fx2 = fx1;

			x1 = a + ratio_a * (b - a);
			fx1 = func(x1);
		}
		else {
			a = x1;

			x1 = x2;
			fx1 = fx2;

			x2 = b + ratio_b * (b - a);
			fx2 = func(x2);
		}

		fprintf(fp, "%u\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\n",
				iter, a, b, b - a, (prev_b - prev_a) / (b - a), x1, x2, fx1, fx2);

		++iter;
	}

	fprintf(fp, "\neps:\t%.14lf\nxmin:\t%.14lf\ni:\t\t%u\nn:\t\t%u\n",
			eps, (a + b) / 2, iter, goldenRatioCalculationNum(iter));

	closeFile(fp);

	return (MethodResult){(a + b) / 2, iter, goldenRatioCalculationNum(iter)};
}

MethodResult fibonacciMethod(double (*func)(double), double a, double b, double eps, const char* filename) {
	FILE* fp = openFile(filename);
	fprintf(fp, "i\ta\tb\t(b - a)\t(prev_b - prev_a)/(b - a)\tx1\tx2\tf(x1)\tf(x2)\n");

	uint32_t n = 0;	
	for (; binet(n + 2) < fabs(b - a) / eps; ++n);

	const double binet_n2 = binet(n + 2);
	const double section_len = b - a;

	double x1 = a + binet(n) / binet_n2 * section_len;
	double x2 = a + b - x1;

	double fx1 = func(x1);
	double fx2 = func(x2);

	double prev_a = a;
	double prev_b = b;

	fprintf(fp, "%u\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\n",
			1, a, b, b - a, (prev_b - prev_a) / (b - a), x1, x2, fx1, fx2);

	uint32_t iter = 2;
	for (; iter != n + 2; ++iter) {
		if (fx1 < fx2) {
			b = x2;

			x2 = x1;
			fx2 = fx1;

			x1 = a + (binet(n - iter + 1) / binet_n2) * section_len;
			fx1 = func(x1);
		}
		else {
			a = x1;

			x1 = x2;
			fx1 = fx2;

			x2 = a + (binet(n - iter + 2) / binet_n2) * section_len;
			fx2 = func(x2);
		}

		fprintf(fp, "%u\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\t%.14lf\n",
				iter, a, b, b - a, (prev_b - prev_a) / (b - a), x1, x2, fx1, fx2);

		prev_a = a;
		prev_b = b;
	}

	fprintf(fp, "\neps:\t%.14lf\nxmin:\t%.14lf\ni:\t\t%u\nn:\t\t%u\n",
			eps, (a + b) / 2, iter, fibonacciCalculationNum(iter));

	closeFile(fp);

	return (MethodResult){(a + b) / 2, iter - 1, fibonacciCalculationNum(iter)};
}

IntervalResult findIntervalMin(double (*func)(double), double x0, double delta, const char* filename) {
	FILE* fp = openFile(filename);
	fprintf(fp, "i\txi\tf(xi)\tinterval\n");

	double fx0 = func(x0);

	if (func(x0 - delta) < fx0)
		delta = -delta;

	double prev_x = x0;
	double x1 = prev_x + delta;
	double h = delta;
	
	double fx1 = func(x1);

	uint32_t iter = 0;
	do {
		prev_x = x0;
		x0 = x1;
		fx0 = fx1;

		h *= 2;
		x1 = prev_x + h;
		fx1 = func(x1);

		fprintf(fp, "%u\t%.14lf\t%.14lf\t[%.14lf, %.14lf]\n", iter, x1, fx1, x0, x1);

		++iter;
	} while (fx1 < fx0);

	closeFile(fp);

	return 0 < delta ? (IntervalResult){prev_x, x1} : (IntervalResult){x1, prev_x};
}
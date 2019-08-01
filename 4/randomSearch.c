#include "randomSearch.h"

#include "rosenbrockMethod.h"

#include <time.h>
#include <math.h>

RSResult simpleRandomSearch(double (*func)(double, double), Point a, Point b, double P, double eps) {
	srand((uint32_t)time(NULL));

	double V = (b.x - a.x) * (b.y - a.y);
	double V_eps = eps * eps;
	double P_eps = V_eps / V;

	uint32_t N = (uint32_t)round(log(1 - P) / log(1 - P_eps));

	double func_min = INFINITY;
	Point point_min = {0, 0};

	uint32_t iter = 0;
	while (iter != N) {
		Point point = {randInRange(a.x, b.x), randInRange(a.y, b.y)};
		
		double func_res = func(point.x, point.y);

		if (func_res < func_min) {
			func_min = func_res;
			point_min = (Point){point.x, point.y};
		}

		++iter;
	}

	return (RSResult){func_min, point_min, iter, iter};
}

RSResult algorithm1(double (*func)(double, double), Point a, Point b, uint32_t m, double eps) {
	srand((uint32_t)time(NULL));

	double func_min = INFINITY;
	Point point_min = {0, 0};
	
	uint32_t iter = 0;
	uint32_t calcs = 0;
	uint32_t miss_search_count = 0;
	while (miss_search_count != m) {
		double point[] = {randInRange(a.x, b.x), randInRange(a.y, b.y)};

		RMResult rosenbrock_res = rosenbrockMethod(func, point, eps);
		calcs += rosenbrock_res.calcs_count;

		if (rosenbrock_res.func_min < func_min) {
			func_min = rosenbrock_res.func_min;
			point_min = (Point){rosenbrock_res.x[0], rosenbrock_res.x[1]};

			miss_search_count = 0;
		}
		else 
			++miss_search_count;

		++iter;
	}

	return (RSResult){func_min, point_min, iter, calcs};
}

RSResult algorithm2(double (*func)(double, double), Point a, Point b, uint32_t m, double eps) {
	srand((uint32_t)time(NULL));
	
	Point point_min = {randInRange(a.x, b.x), randInRange(a.y, b.y)};
	double func_min = func(point_min.x, point_min.y);

	uint32_t iter = 0;
	uint32_t calcs = 1;
	uint32_t miss_search_count = 0;
	while (miss_search_count != m) {
		point_min = (Point){randInRange(a.x, b.x), randInRange(a.y, b.y)};
		++calcs;

		if (func(point_min.x, point_min.y) < func_min && miss_search_count != m) {
			double point[] = {point_min.x, point_min.y};

			RMResult rosenbrock_res = rosenbrockMethod(func, point, eps);
			calcs += rosenbrock_res.calcs_count;

			if (rosenbrock_res.func_min < func_min) {
				func_min = rosenbrock_res.func_min;
				point_min = (Point){rosenbrock_res.x[0], rosenbrock_res.x[1]};

				miss_search_count = 0;
			}
			else 
				++miss_search_count;
		}
		else 
			++miss_search_count;
		
		++iter;
	}

	return (RSResult){func_min, point_min, iter, calcs};
}

RSResult algorithm3(double (*func)(double, double), Point a, Point b, uint32_t m, double eps) {
	srand((uint32_t)time(NULL));

	Point point_min = {randInRange(a.x, b.x), randInRange(a.y, b.y)};
	double func_min = func(point_min.x, point_min.y);

	double point[] = {point_min.x, point_min.y};
	RMResult rosenbrock_res1 = rosenbrockMethod(func, point, eps);

	uint32_t iter = 0;
	uint32_t calcs = 1 + rosenbrock_res1.calcs_count;
	uint32_t miss_search_count = 0;
	while (miss_search_count != m) {
		Point point_step = {randInRange(a.x, b.x), randInRange(a.x, b.x)};

		Point point_temp = {0, 0};
		double func_temp = 0;
		
		do {
			point_temp.x += point_step.x;
			point_temp.y += point_step.y;

			func_temp = func(rosenbrock_res1.x[0] + point_temp.x, rosenbrock_res1.x[1] + point_temp.y);
			++calcs;

			++miss_search_count;
		} while (rosenbrock_res1.func_min < func_temp &&
				 isInRange(a.x, rosenbrock_res1.x[0] + point_temp.x, b.x) &&
				 isInRange(a.y, rosenbrock_res1.x[1] + point_temp.y, b.y) &&
				 miss_search_count != m);

		if (miss_search_count != m) {
			point[0] = rosenbrock_res1.x[0] + point_temp.x;
			point[1] = rosenbrock_res1.x[1] + point_temp.y;

			RMResult rosenbrock_res2 = rosenbrockMethod(func, point, eps);
			calcs += rosenbrock_res2.calcs_count;

			if (rosenbrock_res2.func_min < rosenbrock_res1.func_min) {
				rosenbrock_res1.x[0] = rosenbrock_res2.x[0];
				rosenbrock_res1.x[1] = rosenbrock_res2.x[1];
				rosenbrock_res1.func_min = rosenbrock_res2.func_min;

				func_min = rosenbrock_res1.func_min;
				point_min = (Point){rosenbrock_res1.x[0], rosenbrock_res1.x[1]};

				miss_search_count = 0;
			}
		}

		++iter;
	}

	return (RSResult){func_min, point_min, iter, calcs};
}

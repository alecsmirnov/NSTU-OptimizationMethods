#include "plotter.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define GRID_FILE  "grid"
#define STEPS_FILE "steps.txt"

#define COORDS_COUNT 3

#define X_COORD 0
#define Y_COORD 1

#define STEP_SIZE 0.05

typedef struct ArangeResult {
	double* sequence;

	size_t size;
} ArangeResult;

typedef struct MeshgridResult {
	double** xg;
	double** yg;

	size_t x_size;
	size_t y_size;
} MeshgridResult;

typedef struct PlotSize {
	double x_max;
	double x_min;

	double y_max;
	double y_min;
} PlotSize;

static ArangeResult arange(double begin, double end, double step) {
	ArangeResult result;

	result.size = (size_t)round((end - begin) / step);
	result.sequence = (double*)malloc(sizeof(double) * result.size);

	for (size_t i = 0; i != result.size; ++i)
		result.sequence[i] = begin + step * i;

	return result;
}

static MeshgridResult meshgrid(double* x, double* y, size_t x_size, size_t y_size) {
	MeshgridResult result = (MeshgridResult){NULL, NULL, y_size, x_size};

	result.xg = (double**)malloc(sizeof(double*) * result.x_size);
	result.yg = (double**)malloc(sizeof(double*) * result.x_size);

	for (size_t i = 0; i != result.x_size; ++i)
		result.xg[i] = (double*)malloc(sizeof(double) * result.y_size);

	for (size_t i = 0; i != result.x_size; ++i)
		result.yg[i] = (double*)malloc(sizeof(double) * result.y_size);

	for (size_t i = 0; i != result.x_size; ++i)
		for (size_t j = 0; j != result.y_size; ++j)
			result.xg[i][j] = x[j];

	for (size_t i = 0; i != result.x_size; ++i)
		for (size_t j = 0; j != result.y_size; ++j)
			result.yg[i][j] = y[i];

	return result;
}

static PlotSize findPlotSize() {
	FILE* fp = fopen(TEMP_FILES STEPS_FILE, "r");

	PlotSize plot_size = (PlotSize){DBL_MIN, DBL_MAX, DBL_MIN, DBL_MAX};

	double x, y;
	const uint8_t ELEM_READ_COUNT = 2;
	while (fscanf(fp, "%lf%lf", &x, &y) == ELEM_READ_COUNT) {
		if (plot_size.x_max < x) 
			plot_size.x_max = x;
		if (x < plot_size.x_min) 
			plot_size.x_min = x;

		if (plot_size.y_max < y) 
			plot_size.y_max = y;
		if (y < plot_size.y_min) 
			plot_size.y_min = y;
	}

	plot_size.x_max = round(plot_size.x_max) + 1;
	plot_size.x_min = floor(plot_size.x_min) - 1;

	plot_size.y_max = round(plot_size.y_max) + 1;
	plot_size.y_min = floor(plot_size.y_min) - 1;

	fclose(fp);

	return plot_size;
}

void plotterAddStep(double x, double y) {
	FILE* fp = fopen(TEMP_FILES STEPS_FILE, "a");

	fprintf(fp, "%lf\t%lf\n", x, y);

	fclose(fp);
}

void plotterClearSteps() {
	fclose(fopen(TEMP_FILES STEPS_FILE, "w"));
}

static void makeGrid(double (*func)(double, double)) {
	PlotSize plot_size = findPlotSize();

	ArangeResult x = arange(plot_size.x_min, plot_size.x_max, STEP_SIZE);
	ArangeResult y = arange(plot_size.y_min, plot_size.y_max, STEP_SIZE);

	MeshgridResult mg = meshgrid(x.sequence, y.sequence, x.size, y.size);
	double** zg = (double**)malloc(sizeof(double*) * mg.x_size);
	for (size_t i = 0; i != mg.x_size; ++i)
		zg[i] = (double*)malloc(sizeof(double) * mg.y_size);

	for (size_t i = 0; i != mg.x_size; ++i)
		for (size_t j = 0; j != mg.y_size; ++j)
			zg[i][j] = func(mg.xg[i][j], mg.yg[i][j]);

	char coords_alias[COORDS_COUNT] = {'x', 'y', 'z'};
	double** grids[COORDS_COUNT] = {mg.xg, mg.yg, zg};
	size_t sizes[COORDS_COUNT][2] = {{mg.x_size, mg.y_size},{mg.x_size, mg.y_size},{mg.x_size, mg.y_size}};

	for (uint8_t i = 0; i != COORDS_COUNT; ++i) {
		char grid_filename[FILENAME_MAX];
		sprintf(grid_filename, TEMP_FILES GRID_FILE "_%c.txt", coords_alias[i]);

		FILE* fp = fopen(grid_filename, "w");

		for (size_t j = 0; j != sizes[i][X_COORD]; ++j) {
			for (size_t k = 0; k != sizes[i][Y_COORD]; ++k)
				fprintf(fp, "%lf ", grids[i][j][k]);

			fprintf(fp, "\n");
		}

		fclose(fp);
	}

	free(x.sequence);
	free(y.sequence);

	for (uint8_t i = 0; i != COORDS_COUNT; ++i) {
		for (size_t j = 0; j != sizes[i][X_COORD]; ++j)
			free(grids[i][j]);

		free(grids[i]);
	}
}

void plotterMakePicture(double (*func)(double, double), const char* funcname, const char* result_path) {
	makeGrid(func);

	char* python_str = (char*)malloc(sizeof(char) *
					   (strlen("python plotter.py ") + strlen(TEMP_FILES) + strlen(funcname) + strlen(result_path) + 4));
	sprintf(python_str, "python plotter.py %s %s %s\0", TEMP_FILES, funcname, result_path);

	system(python_str);

	free(python_str);
}

void plotterClearData() {
	clearSteps();

	fclose(fopen(TEMP_FILES GRID_FILE "_x.txt", "w"));
	fclose(fopen(TEMP_FILES GRID_FILE "_y.txt", "w"));
	fclose(fopen(TEMP_FILES GRID_FILE "_z.txt", "w"));
}
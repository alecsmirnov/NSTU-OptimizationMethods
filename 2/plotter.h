#ifndef PLOTTER_H
#define PLOTTER_H

#define TEMP_FILES "plotter_temp/"

void plotterAddStep(double x, double y);
void plotterClearSteps();
void plotterMakePicture(double (*func)(double, double), const char* funcname, const char* result_path);
void plotterClearData();

#endif
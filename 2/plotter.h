#ifndef PLOTTER_H
#define PLOTTER_H

#define TEMP_FILES "plotter_temp/"

void addStep(double x, double y);
void clearSteps();
void makePicture(double (*func)(double, double), const char* funcname, const char* result_path);
void clearData();

#endif
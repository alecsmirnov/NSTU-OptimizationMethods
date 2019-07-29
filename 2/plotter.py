import matplotlib.pyplot as plt
import sys
import re

TEMP_DIR = 'plotter_temp/'

ARRAY_FIRST =  0
ARRAY_LAST  = -1

LINE_WIDTH 		 = 1
LINE_COLOR 		 = 'blue'
LINE_BEGIN_COLOR = 'orange'
LINE_END_COLOR 	 = 'green'

FUNCNAME_LABEL_FONT_SIZE = 18
LABEL_FONT_SIZE 		 = 12

LEVELS_COUNT = 26

def getCoordSteps():
	x_steps = []
	y_steps = []
	with open(TEMP_DIR + 'steps.txt') as file:
 		for line in file:
 			x, y = line.split()
 			x_steps.append(float(x))
 			y_steps.append(float(y))
	return x_steps, y_steps

def getGrids():
	xg = list(map(str.split, open(TEMP_DIR + 'grid_x.txt')))
	yg = list(map(str.split, open(TEMP_DIR + 'grid_y.txt')))
	zg = list(map(str.split, open(TEMP_DIR + 'grid_z.txt')))
	return xg, yg, zg

def setPlotSettings():
	plt.title(funcname, fontsize = FUNCNAME_LABEL_FONT_SIZE)
	plt.xlabel("X", fontsize = LABEL_FONT_SIZE)
	plt.ylabel("Y", fontsize = LABEL_FONT_SIZE, rotation = 0)

def drawLine(x_steps, y_steps):
	plt.plot(x_steps, y_steps, linewidth = LINE_WIDTH, color = LINE_COLOR)
	plt.scatter(x_steps[ARRAY_FIRST], y_steps[ARRAY_FIRST], color = LINE_BEGIN_COLOR)
	plt.scatter(x_steps[ARRAY_LAST], y_steps[ARRAY_LAST], color = LINE_END_COLOR)

def makePicture(funcname, result_path):
	setPlotSettings()
	x_steps, y_steps = getCoordSteps()
	drawLine(x_steps, y_steps)
	xg, yg, zg = getGrids()
	plt.contour(xg, yg, zg, LEVELS_COUNT)
	plt.savefig(result_path + funcname + ".png")

if __name__ == "__main__":	
	TEMP_DIR = sys.argv[1]
	funcname = re.sub('\.txt$', '', sys.argv[2])
	result_path = sys.argv[3]
	makePicture(funcname, result_path)
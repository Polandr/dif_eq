#!/usr/bin/python

import matplotlib.pyplot as plt
from math import sqrt
import sys
import time

plt.ion()
plt.show()

y_max = 0
y_min = 0

x_axis = sys.stdin.readline()
x_values = [float(x) for x in x_axis.split()]
#print x_values

itr = 0

while 1:
	try:
		line = sys.stdin.readline()
		if not line:
			break

		T = float(line)

		line_re = sys.stdin.readline()
		if not line_re:
			break
		line_im = sys.stdin.readline()
		if not line_im:
			break

		y_values_re = [float(x) for x in line_re.split()]
		y_values_im = [float(x) for x in line_im.split()]
		y_values_abs = [sqrt(y_values_re[i]*y_values_re[i] + y_values_im[i]*y_values_im[i]) for i in range(0,len(y_values_re))]

		if (itr == 0 or y_max < max(y_values_re) or y_max < max(y_values_im) or y_max < max(y_values_abs)):
			max_re = max(y_values_re)
			max_im = max(y_values_im)
			max_abs = max(y_values_abs)
			y_max = max(max_re,max_im)
			y_max = max(y_max,max_abs)


		if (itr == 0 or y_min > min(y_values_re) or y_min > min(y_values_im) or y_min > min(y_values_abs)):
			min_re = min(y_values_re)
			min_im = min(y_values_im)
			min_abs = min(y_values_abs)
			y_min = min(min_re,min_im)
			y_min = min(y_min,min_abs)


		plt.clf()
		plt.axis([x_values[0], x_values[-1], y_min, y_max])
		plt.plot(x_values, y_values_re, 'r', antialiased = True)
		plt.plot(x_values, y_values_im, 'g', antialiased = True)
		plt.plot(x_values, y_values_abs, 'b', antialiased = True)
		print "visualize [T = %f]" % T
		plt.draw()

		time.sleep( 0.00001 )

		itr += 1

	except KeyboardInterrupt:
		sys.exit()

while 1:
	try:
		sys.stdin.readline()
	except KeyboardInterrupt:
		sys.exit()

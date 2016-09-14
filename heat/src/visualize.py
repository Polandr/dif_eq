#!/usr/bin/python

import matplotlib.pyplot as plt
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

		line = sys.stdin.readline()
		if not line:
			break

		y_values = [float(x) for x in line.split()]
		#print y_values

		if (itr == 0 or y_max < max(y_values)):
			y_max = max(y_values)
		if (itr == 0 or y_min > min(y_values)):
			y_min = min(y_values)

		#plt.clf()
		plt.axis([x_values[0], x_values[-1], y_min, y_max])
		plt.plot(x_values, y_values, antialiased = True)
		print "visualize [T = %.4f]" % T
		plt.draw()

		time.sleep( 0.1 )

		itr += 1

	except KeyboardInterrupt:
		sys.exit()

while 1:
	try:
		sys.stdin.readline()
	except KeyboardInterrupt:
		sys.exit()

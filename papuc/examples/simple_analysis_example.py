
from matplotlib import pyplot as plt
import numpy as np


#####################
from papuc import isoluminant_uniform_spline_colormap

#### define a color map with knot points
L_max = 85 #   0   1    2    3    4    5
a_knots =  [  -6,  16,  15,  0,  -22, -22]
b_knots =  [ -15, -13,  -3,  0,   -1, -10]
my_map = isoluminant_uniform_spline_colormap(
	a_knots, b_knots, L_max,
	maximize_radius = False) # this will expand the shape to the bounds of sRGB1 color space



###### now for a basic path plot
from papuc.analysis import cyclic_colormap_test, plot_knots_on_isoluminant_slice

fig, ax = plt.subplots()
plot_knots_on_isoluminant_slice(ax, my_map)
fig.savefig('colormap_path_knots.png')

# the size can be maximized later too!
fig, ax = plt.subplots()
my_map.maximize_radius()
plot_knots_on_isoluminant_slice(ax, my_map)
fig.savefig('colormap_path_knots_maximized.png')



# Enable this if you want to do a more rigorous analysis, it make take a minute to run though
cyclic_colormap_test(my_map)
fig.savefig('colormap_analyzed.png')

plt.show()

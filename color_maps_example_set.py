

from colorspacious import cspace_convert
from colorsys import hsv_to_rgb
from analysis import cyclic_colormap_test
import numpy as np

jab_color = cspace_convert(hsv_to_rgb(0.0,1.0,1.0),"sRGB1",  "CAM02-UCS")
HSV_ab_angle0 = np.arctan2(jab_color[2], jab_color[1])
print(jab_color, HSV_ab_angle0 * 180/np.pi)


test_all = False





name = 'bean5'
L_max = 75 #   0   1    2    3    4    5
a_knots =  [   0,  23,  26,  0, -23, -22]
b_knots =  [ -24, -19,  -3,  0,  25, -15]
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False)

cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = True, ab_grid =512)




## for CAM02-UCS
L_max = 73
radius = 26
name = 'circle-%.1f-%.1f' %(radius,L_max)
center = (0,0)
theta = np.arange(0, 2*np.pi, 2*np.pi/8)
a_knots = center[0] + radius*np.cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*np.sin(theta + HSV_ab_angle0)
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False)

#cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = True, ab_grid =512)









## literaly random color cycle
name = 'random'
L_max = 73
from numpy.random import randint
knots = 7
# good for L*=75
amin, amax = -30, 30
bmin, bmax = -30, 30
a_knots = randint(amin, amax, size = knots)
b_knots = randint(bmin, bmax, size = knots)
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name)

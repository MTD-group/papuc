

from colorspacious import cspace_convert
from colorsys import hsv_to_rgb
from analysis import cyclic_colormap_test
from numpy import *

lab_color = cspace_convert(hsv_to_rgb(0.0,1.0,1.0),"sRGB1",  "CIELab")
HSV_ab_angle0 = arctan2(lab_color[2], lab_color[1])
print(lab_color, HSV_ab_angle0 * 180/pi)


test_all = False




name = 'bean1'
L_max = 75
a_knots =     [ -30, 40, 30,  -11, -60]
b_knots =  [  -20, -25 ,  -10,  9, 60]
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, color_space = "CIELab")


name = 'bean2'
L_max = 75
a_knots =     [ -30, 45, 30,  0, -60]
b_knots =  [  -20, -30 ,  -10, 0, 60]
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, color_space = "CIELab")



name = 'bean3'
L_max = 75
a_knots =     [ -30, 30, 40, -15,   -60]
b_knots =  [  -20, -35 , -10, 15,  40]
if test_all: cyclic_colormap_test(a_knots, b_knots, L_max, name, color_space = "CIELab")



name = 'bean4'
L_max = 75 #   0   1    2    3    4    5
a_knots =  [ -13,  48,  35,  0, -65, -33]
b_knots =  [ -36, -33,  -5,  0,  50, -20]
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False, color_space = "CIELab")

#cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = True)


###### perimeter
name = 'perimeter1'
L_max = 74
#                    0   1,   2,    3,   4,   5,    6,   7,  8,    9
a_knots = array( [  34,  27,  12, -50, -68, -50,  -26,  -7,  42,  50])
b_knots = array( [  29,  62,  75,  70,  52,   5,  -28, -39, -37, -20])
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False, color_space = "CIELab")

#cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False)


### a circle


## for Lab
L_max = 74
radius = 40
name = 'circle-%.1f-%.1f' %(radius,L_max)
center = (0,0)
theta = arange(0, 2*pi, 2*pi/8)
a_knots = center[0] + radius*cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*sin(theta + HSV_ab_angle0)
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False, color_space = "CIELab")

cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = True, ab_grid =512, color_space = "CIELab")






### a darker circle
L_max = 50
radius = 29
name = 'circle-%.1f-%.1f' %(radius,L_max)
center = (0,0)
theta = arange(0, 2*pi, 2*pi/8)
a_knots = center[0] + radius*cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*sin(theta + HSV_ab_angle0)
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False, color_space = "CIELab")



### a lighter circle
L_max = 85
radius = 21
name = 'circle-%.1f-%.1f' %(radius,L_max)
center = (0,0)
theta = arange(0, 2*pi, 2*pi/8)
a_knots = center[0] + radius*cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*sin(theta + HSV_ab_angle0)
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False, color_space = "CIELab")

#cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = True)


## literaly random color cycle
name = 'random'
L_max = 75
from numpy.random import randint
knots = 10
# good for L*=75
amin, amax = -80, 65
bmin, bmax = -45, 85
a_knots = randint(amin, amax, size = knots)
b_knots = randint(bmin, bmax, size = knots)
if test_all:cyclic_colormap_test(a_knots, b_knots, L_max, name, color_space = "CIELab")

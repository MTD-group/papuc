
from colorspacious import cspace_convert
from colorsys import hsv_to_rgb
from papuc.core import isoluminant_uniform_spline_colormap
import numpy as np


colormaps = {}


lab_color = cspace_convert(hsv_to_rgb(0.0,1.0,1.0),"sRGB1",  "CIELab")
HSV_ab_angle0 = np.arctan2(lab_color[2], lab_color[1])
#print(lab_color, HSV_ab_angle0 * 180/np.pi)







name = 'bean1'
L_max = 75
a_knots =     [ -30, 40, 30,  -11, -60]
b_knots =  [  -20, -25 ,  -10,  9, 60]
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map



name = 'bean2'
L_max = 75
a_knots =     [ -30, 45, 30,  0, -60]
b_knots =  [  -20, -30 ,  -10, 0, 60]
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map



name = 'bean3'
L_max = 75
a_knots =     [ -30, 30, 40, -15,   -60]
b_knots =  [  -20, -35 , -10, 15,  40]
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map



name = 'bean4'
L_max = 75 #   0   1    2    3    4    5
a_knots =  [ -13,  48,  35,  0, -65, -33]
b_knots =  [ -36, -33,  -5,  0,  50, -20]
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map




###### perimeter
name = 'perimeter1'
L_max = 74
#                       0   1,   2,    3,   4,   5,    6,   7,  8,    9
a_knots = np.array( [  34,  27,  12, -50, -68, -50,  -26,  -7,  42,  50])
b_knots = np.array( [  29,  62,  75,  70,  52,   5,  -28, -39, -37, -20])
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map




### a circle


## widest for Lab
L_max = 74
radius = 40
name = 'circle-%.1f-%.1f' %(radius,L_max)
center = (0,0)
theta = np.arange(0, 2*np.pi, 2*np.pi/8)
a_knots = center[0] + radius*np.cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*np.sin(theta + HSV_ab_angle0)
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map
colormaps['default'] = map





### a darker circle
L_max = 50
radius = 29
name = 'circle-%.1f-%.1f' %(radius,L_max)
center = (0,0)
theta = np.arange(0, 2*np.pi, 2*np.pi/8)
a_knots = center[0] + radius*np.cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*np.sin(theta + HSV_ab_angle0)
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map



### a lighter circle
L_max = 85
radius = 21
name = 'circle-%.1f-%.1f' %(radius,L_max)
center = (0,0)
theta = np.arange(0, 2*np.pi, 2*np.pi/8)
a_knots = center[0] + radius*np.cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*np.sin(theta + HSV_ab_angle0)
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map




## literaly random color cycle
name = 'random'
L_max = 75
from numpy.random import randint
knots = 7
# good for L*=75
amin, amax = -80, 65
bmin, bmax = -45, 85
a_knots = randint(amin, amax, size = knots)
b_knots = randint(bmin, bmax, size = knots)
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = "CIELab")
colormaps[name] = map

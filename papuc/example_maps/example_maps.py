
from colorspacious import cspace_convert
from colorsys import hsv_to_rgb
from papuc.core import isoluminant_uniform_spline_colormap
import numpy as np


colormaps = {}






## for CAM02-UCS
L_max = 73
radius = 26
name = 'circle-%.1f-%.1f' %(radius,L_max)
jab_color = cspace_convert(hsv_to_rgb(0.0,1.0,1.0),"sRGB1",  "CAM02-UCS")
HSV_ab_angle0 = np.arctan2(jab_color[2], jab_color[1])
#print(jab_color, HSV_ab_angle0 * 180/np.pi)
center = (0,0)
theta = np.arange(0, 2*np.pi, 2*np.pi/8)
a_knots = center[0] + radius*np.cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*np.sin(theta + HSV_ab_angle0)
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max)
colormaps[name] = map
colormaps['default'] = map





L_max = 73
radius = 25 # will get rescaled later
name = 'max-circle-%i' %(L_max)
jab_color = cspace_convert(hsv_to_rgb(0.0,1.0,1.0),"sRGB1",  "CAM02-UCS")
HSV_ab_angle0 = np.arctan2(jab_color[2], jab_color[1])
#print(jab_color, HSV_ab_angle0 * 180/np.pi)
center = (0,0)
theta = np.arange(0, 2*np.pi, 2*np.pi/8)
a_knots = center[0] + radius*np.cos(theta + HSV_ab_angle0)
b_knots = center[1] + radius*np.sin(theta + HSV_ab_angle0)
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, maximize_radius = True, verbose = False)
colormaps[name] = map
if __name__ == '__main__':
    print('optimal radius found', np.sqrt(np.mean(map.a_knots**2 + map.b_knots**2)))






#### an attempt at avoiding the red-green spectrum for r-g colorblindness
name = 'bean5'
L_max = 75 #   0   1    2    3    4    5
a_knots =  [   0,  23,  26,  0, -23, -22]
b_knots =  [ -24, -19,  -3,  0,  25, -15]
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, maximize_radius = True)
colormaps[name] = map




#### literaly random color cycle
name = 'random'
L_max = 73
from numpy.random import randint
knots = 7
# good for L*=75
amin, amax = -30, 30
bmin, bmax = -30, 30
a_knots = randint(amin, amax, size = knots)
b_knots = randint(bmin, bmax, size = knots)
map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max)
colormaps[name] = map


from numpy import arange, meshgrid, ones,zeros, linspace, pi
from numpy import any as np_any

num_theta_levels = 5
num_z_levels = 8


dtheta = 1.0/num_theta_levels
theta_levels = arange(dtheta/2.0, 1.0, dtheta)
dz = 1.0/num_z_levels
z_levels = arange(dz/2.0, 1.0, dz)


#sequence = my_map(theta_levels, ones(theta_levels.shape))
#deut_sequence = cspace_convert(sequence, cvd_space, "sRGB1")

#theta_grid_2D, z_grid_2D =  meshgrid(theta_levels, z_levels, sparse=False, indexing='ij')
#image = my_map(theta_grid_2D, z_grid_2D)
#deut_image = cspace_convert(image, cvd_space, "sRGB1")


from tools import isoluminant_uniform_circle_colormap

my_map = isoluminant_uniform_circle_colormap(radius = 80, L_max = 75)

L_max_values = arange(50,80,1)
r_max_values = []
for L_max in L_max_values:
	my_map.L_max = L_max
	my_map.maximize_radius(cut_off = 0.1)
	r_max_values.append(my_map.radius)
	print (my_map.L_max, my_map.radius)





from matplotlib.pylab import *

plot(L_max_values, r_max_values)

show()


from matplotlib.pyplot import figure, quiver, show, imshow, imsave, subplots
from numpy import meshgrid, sin, cos, linspace, savetxt, loadtxt, pi, sqrt, angle, arctan2



######### this first part synthesizes some data to try
npts = 16

xmax = 2
ymax = 1

X, Y = meshgrid(linspace(0, xmax , xmax * npts), linspace(0, ymax , ymax*npts))
U = sin(2*pi*Y)
V = cos(2*pi*X)


savetxt('x_vectors.txt',U)
savetxt('y_vectors.txt',V)

figure('Simple Analytical Test')
quiver(X, Y, U, V, units='width')





###################### read the data

U = loadtxt('x_vectors.txt')
V = loadtxt('y_vectors.txt')


# you could also use our example data
# z_array = loadtxt('example_data/magnitudes.txt')
# theta_array = loadtxt('example_data/theta_angles.txt')


angle = arctan2(V,U)
magnitude = sqrt(V**2 + U**2)

magnitude_norm = magnitude/magnitude.max()
angle_norm     = angle/(2*pi)



from papuc import isoluminant_uniform_spline_colormap
## in CAM02-UCS coordinates be default
L_max = 73
radius = 26
theta0 = pi/8.0 # phase shift the colormap
theta = linspace(theta0 , theta0 + 2*pi, 8, endpoint=False)

my_map = isoluminant_uniform_spline_colormap(  # initialize the colormap
	a_knots = radius*cos(theta),
	b_knots = radius*sin(theta),
	L_max = L_max)
my_map.name = 'circle-%.1f-%.1f' %(radius,L_max)

image = my_map(angle_norm, magnitude_norm)
imsave('simple_test_image.png', image, origin = 'lower')

### that's all you need!
# this is what it looks like
figure('From File')
imshow(image, origin = 'lower',      extent = (0, xmax, 0, ymax) , interpolation ='None' )
quiver(X, Y, U, V, units='width')





### this is for looking at your colormap's color wheel and color path
#figure('Color Wheel and Color Path')
fig, axes = subplots(ncols=2, num='Color Wheel and Color Path')
from papuc import plot_knots_on_isoluminant_slice, plot_colorwheel
plot_colorwheel(axes[0], my_map)
plot_knots_on_isoluminant_slice(axes[1], my_map)
fig.tight_layout(pad= 0.1)

show()

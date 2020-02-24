from matplotlib import pyplot as plt
import numpy as np
##### Synthesizes some data to plot
npts = 16
xmax = 2
ymax = 1
X, Y = np.meshgrid(
    np.linspace(0, xmax , xmax * npts),
    np.linspace(0, ymax , ymax * npts))
U = np.sin(2*np.pi*Y)
V = np.cos(2*np.pi*X)
angle = np.arctan2(V,U)
magnitude = np.sqrt(V**2 + U**2)
magnitude_norm = magnitude/magnitude.max()

###### Now for visualization 
from papuc.example_maps import colormaps
my_map = colormaps['default']

## make an sRGB1 image from the data
image = my_map(angle, magnitude_norm)

## this is what it looks like with a quiver plot for comparison
plt.imshow(image, origin = 'lower',  extent = (0, xmax, 0, ymax) )
plt.quiver(X, Y, U, V, units='width')

## save the image
plt.imsave('test_image.png', image, origin = 'lower')

## save the combined figure
plt.savefig('test_figure.png')

## this is for looking at your colormap's color wheel
fig, ax = plt.subplots()
from papuc.analysis import plot_colorwheel, UCS_cone_3D
plot_colorwheel(ax, my_map)
fig.tight_layout(pad= 0.1)
plt.savefig('test_colorwheel.png')

## a nice 3D view
fig = plt.figure()
ax = fig.gca(projection='3d')
UCS_cone_3D(ax,  my_map)
fig.savefig('test_colormap_3D_view.png')

plt.show()

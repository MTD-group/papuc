# Phase-Amplitude Perceptually Uniform Colormap Designer (papuc) ##
Papuc is a small library designed for creating perceptually uniform color mappings of vector fields where vector orientations/phase angles are mapped onto hue and vector magnitudes/amplitudes are mapped onto lightness. papuc is designed to integrate into the Matplotlib/NumPy ecosystem and relies on [colorspacious](https://colorspacious.readthedocs.io/en/latest/) for colorspace transformations. 

## Why? ##
Simply put, many people use the [HSV color wheel](https://en.wikipedia.org/wiki/HSL_and_HSV) for colormaping phase and amplitude. The HSV color wheel is handy for artists picking color palletes, but has serious failings as a color mapping for data visualization.   

## Citing Papuc ##
Our paper exploring the concepts, the improvements over HSV, and a usecase are [here](Will_be_on_ArXiv_soon)

## Basic Usage ##
If you just want to use a resonably good color map on your angle and magnitude data, our usage_example.py should be all you need. 
```
>>> from numpy import pi
>>> from papuc.example_maps import colormaps
>>> my_map = colormaps['default']


## A Feature Complete Example ##
We start by creating some test data:
```python
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
```
Now we make an sRGB1 image (r, g, and b scaled from 0.0 to 1.0) for our data:
```python
###### Now for visualization 
from papuc.example_maps import colormaps
my_map = colormaps['default']

## make an sRGB1 image from the data
image = my_map(angle, magnitude_norm)
```
This next part will plot the colormapped data with a quiver plot of the same vector field and show the colorwheel. 
```python
plt.imshow(image, origin = 'lower',  extent = (0, xmax, 0, ymax) )
plt.quiver(X, Y, U, V, units='width')

## save the image
plt.imsave('test_image.png', image, origin = 'lower')

## save the combined figure
plt.savefig('test_figure.png')

## this is for looking at your colormap's color wheel
fig, ax = plt.subplots()
from papuc.analysis import plot_colorwheel
plot_colorwheel(ax, my_map)
fig.tight_layout(pad= 0.1)
plt.savefig('test_colorwheel.png')

plt.show()
```
Which looks like this:

![test_figure](/papuc/examples/test_figure.png)

![test_colorwheel](/papuc/examples/test_colorwheel.png)

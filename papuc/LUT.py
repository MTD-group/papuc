
import matplotlib.pyplot as plt
import numpy as np


def create_LUT(my_map, npts = 128, save_to_file = True):

	xmax = 1
	ymax = 2*np.pi

	magnitude_norm, angle = np.meshgrid(np.linspace(0, xmax , xmax * npts), np.linspace(0, ymax , ymax*npts))

	image = my_map(angle, magnitude_norm)
	if save_to_file:
		plt.imsave(my_map.name+'_LUT_image.png', image, origin = 'lower')
		np.savetxt(my_map.name+'_LUT_red.txt',  image[:,:,0])
		np.savetxt(my_map.name+'_LUT_green.txt',image[:,:,1])
		np.savetxt(my_map.name+'_LUT_blue.txt', image[:,:,2])

	return image



if __name__ == '__main__':

	name = 'default'
	from papuc.example_maps import colormaps

	my_map = colormaps[name]
	image = create_LUT(my_map, npts = 32, save_to_file = False)

	plt.imshow(image, origin = 'lower',      extent = (0, 1.0, 0, np.pi*2) , interpolation ='None' )
	plt.show()

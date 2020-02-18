
import matplotlib.pyplot as plt
import numpy as np
from papuc.LUT import create_LUT



name = 'default'
from papuc.example_maps import colormaps

my_map = colormaps[name]
image = create_LUT(my_map, npts = 32, save_to_file = False)

plt.imshow(image, origin = 'lower',      extent = (0, 1.0, 0, np.pi*2) , interpolation ='None' )
plt.show()

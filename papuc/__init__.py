from .core import isoluminant_uniform_spline_colormap,  isoluminant_uniform_circle_colormap
from .analysis import plot_knots_on_isoluminant_slice, plot_colorwheel
from .example_maps import colormaps



from .core import isoluminant_uniform_spline_colormap, periodic_spline, periodic_spline_test, default_colorspace
from .analysis import plot_knots_on_isoluminant_slice, plot_colorwheel, cyclic_colormap_test


# these functions are exposed to the user for development purposes, you don't need them for using a prexisting color map
from .analysis import (plot_analytic_vortex_test_data,
	plot_small_wave_magnitude_test,
	plot_small_wave_angle_test,
	plot_perceptual_derivative,
	plot_arc_length_vs_spline_parameter_t)
	

from .plot_3d     import UCS_cone_3D
from .plot_hsv_3d import hsv_UCS_cone_3D

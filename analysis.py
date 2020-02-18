

#######################
#from numpy import array, pi, linspace, mod,  append, arange, sqrt, diff, mean, clip, zeros, any, all, cos, sin, floor, meshgrid, ones, arctan2
import numpy as np
from scipy.interpolate import splrep, splev

from colorspacious import cspace_convert

from papuc.tools import periodic_spline, periodic_spline_test, isoluminant_uniform_spline_colormap
from papuc.tools import HSV_colormap, default_colorspace

if False:
	periodic_spline_test()


def plot_knots_on_isoluminant_slice(ax, my_map, ab_grid=512):
	na_grid = nb_grid = ab_grid
	amin, amax = -50, 50
	bmin, bmax = -50, 50

	if my_map.color_space == "CIELab":
			amin, amax = -128, 128
			bmin, bmax = -128, 128

	color_grid = True

	ax.set_title( my_map.color_space)


	### plot the smooth spline path
	if my_map.min_theta_pts != None:
		theta_plot = np.linspace(0,2*np.pi, my_map.min_theta_pts*4)
		ax.plot(
				my_map.a_func_of_theta(theta_plot),
				my_map.b_func_of_theta(theta_plot), color ='k')

		### plot the equally spaced points from the minimum v_refinement
		theta_max_space = np.linspace(0, 2*np.pi, my_map.min_theta_pts)
		ax.plot(
				my_map.a_func_of_theta(theta_max_space),
				my_map.b_func_of_theta(theta_max_space), marker = 'x', color = 'k',linestyle = '' )
		### plot the spline knots
		ax.plot(my_map.a_knots,my_map.b_knots, marker = 'o', linestyle = '', mfc = 'white', mec = 'black', markersize = 10)
		for i in range(len(my_map.a_knots)):
			ax.text(
						my_map.a_knots[i],
						my_map.b_knots[i], '%i'%i, ha = 'center', va = 'center', fontsize = 8)
	else:
		theta_plot = np.linspace(0,2*np.pi, 360)
		ax.plot(
			my_map.a_func_of_theta(theta_plot),
			my_map.b_func_of_theta(theta_plot), color ='k')

		ax.plot(my_map.a_func_of_theta(my_map.theta0),my_map.b_func_of_theta(my_map.theta0), marker = 'o', linestyle = '', mfc = 'white', mec = 'black', markersize = 10)
		i = 0
		ax.text(
					my_map.a_func_of_theta(my_map.theta0),
					my_map.b_func_of_theta(my_map.theta0), '%i'%i, ha = 'center', va = 'center', fontsize = 8)

	ax.set_xlabel('a')
	ax.set_ylabel('b')
	ax.axis('equal')

	if color_grid:

		use_image = True
		if use_image:
			print('making image')
			da = (amax-amin)/na_grid
			db = (bmax-bmin)/nb_grid
			a_points = np.linspace(amin+da/2, amax-da/2, na_grid)
			b_points = np.linspace(bmin+db/2, bmax-db/2, nb_grid)

			a_grid, b_grid = np.meshgrid(a_points,b_points, indexing = 'xy')

			lab_map = np.ones((na_grid, nb_grid, 3))*my_map.L_max

			lab_map[:,:,1] = a_grid
			lab_map[:,:,2] = b_grid

			bg_color = np.array([1, 1, 1, 0])

			sRGB1_map = np.zeros((na_grid, nb_grid, 4)) # default is transparent
			sRGB1_map[:,:,0:3] = cspace_convert(lab_map, my_map.color_space , "sRGB1")
			for a_index in range(na_grid):
				for b_index in range(nb_grid):
					sRGB1_color = sRGB1_map[a_index,b_index,0:3]
					if  all(0<sRGB1_color) and all(sRGB1_color<1):
						sRGB1_map[a_index,b_index,3] = 1
					else:
						sRGB1_map[a_index,b_index] = bg_color

			ax.imshow(sRGB1_map, extent=(amin,amax,bmin,bmax),
				 origin = 'lower', interpolation = 'bicubic') #'none', 'nearest', 'bilinear', 'bicubic'

		else: # This will use slow patches
			import matplotlib.patches as mpatches
			a_edges, da = np.linspace(amin,amax, na_grid+1, retstep=True)
			b_edges, db = np.linspace(bmin,bmax, nb_grid+1, retstep=True)

			sRGB1_edgecolor = cspace_convert([100-L_max,0,0], my_map.color_space , "sRGB1")

			for a_index in range(na_grid):
				for b_index in range(nb_grid):
					a = a_edges[a_index]+0.5*da
					b = b_edges[b_index]+0.5*db

					lab_color = [L_max, a, b]
					sRGB1_color = cspace_convert(lab_color, my_map.color_space , "sRGB1")

					if  all(0<sRGB1_color) and all(sRGB1_color<1):

						ax.add_patch(mpatches.Rectangle(
								xy = (a_edges[a_index], b_edges[b_index]),
								width = da,
								height = db,
								color = sRGB1_color))#,
								#linestyle = '-',linewidth = 4.0,
								#edgecolor = 'k'))



def plot_arc_length_vs_spline_parameter_t(ax, my_map):

	arc_len_in_t = np.linspace(0.0,1.0, my_map.min_t_pts)
	arc_len = my_map.arc_length * my_map.func_theta_of_t(arc_len_in_t)/(np.pi*2)
	arc_len[-1] = my_map.arc_length
	ax.plot(arc_len_in_t, arc_len )

	theta_max_space = np.linspace(0, 2*np.pi ,my_map.min_theta_pts)
	t_at_theta_max_space = my_map.func_t_of_theta(theta_max_space)
	t_at_theta_max_space[-1] = 1.0
	ax.plot(t_at_theta_max_space, theta_max_space* my_map.arc_length/(np.pi*2), marker = '+', linestyle = '', color = 'k')
	ax.set_xlabel('Spline Parametric, t')
	ax.set_ylabel('Arc Length')
	ax.set_xlim(0,1.0)
	ax.set_ylim(0,ax.get_ylim()[1] )


def plot_perceptual_derivative(ax, my_map, use_delta_instead = False):
	for theta in my_map.theta_knots:
		ax.axvline(theta, color = 'grey')

	dE = my_map.compute_min_theta_points_for_uniformity()
	dv = 2*np.pi/(len(dE))
	v_points_for_center_diff = np.arange(dv/2.0, 2*np.pi, dv)

	if use_delta_instead:
		y =  dE
		ax.set_ylabel(r'$\Delta E$')
	else:
		y = dE/dv
		ax.set_ylabel(r'$\frac{dE}{d \theta}$')

	ax.plot( v_points_for_center_diff, y, marker = 'o', markersize = 2, label = '%i pts'%(len(dE)+1))
	ax.set_xlabel(r'Angle, $\theta$,')
	ax.set_xlim((0,np.pi*2))
	ax.legend(handlelength = 1.0)





def plot_colorwheel(ax, my_map, num_theta_levels = 32, num_z_levels = 16, center = (0,0), radius=1.0, with_side_bar = True, with_deuteranomaly = False):
	import matplotlib.patches as mpatches

	############
	dtheta = 1.0/num_theta_levels
	theta_levels = 2*np.pi* np.arange(dtheta/2.0, 1.0, dtheta)
	dtheta = dtheta * 2*np.pi
	dz = 1.0/num_z_levels
	z_levels = np.arange(dz/2.0, 1.0, dz)
	#################

	cvd_space = {"name": "sRGB1+CVD",
		"cvd_type": "deuteranomaly",
		"severity": 50}

	sRGB1_sequence = my_map(theta_levels, np.ones(theta_levels.shape))
	deut_sequence = np.clip( cspace_convert(sRGB1_sequence, cvd_space, "sRGB1"), 0, 1) # keep it in bounds

	theta_grid_2D, z_grid_2D =  np.meshgrid(theta_levels, z_levels, sparse=False, indexing='ij')
	sRGB1_radial_image = my_map(theta_grid_2D, z_grid_2D, verbose = False)
	deut_radial_image = np.clip( cspace_convert(sRGB1_radial_image, cvd_space, "sRGB1"), 0, 1) # keep it in bounds

	xc = center[0]
	yc = center[1]

	for ti in range(num_theta_levels):

		if with_side_bar:
			if with_deuteranomaly == False:
				color = sRGB1_sequence[ti]
			else:
				color =  deut_sequence[ti]

			ax.add_patch(mpatches.Rectangle(
						xy = (xc +1.1*radius ,  yc + radius*(2.0*theta_levels[ti]/(2*np.pi) - dtheta/(2*np.pi) - 1.0 ) ), # bottom is dtheta/2 but we mulitply by 2 radii
						width = 0.2,
						height = radius*dtheta*2/(2*np.pi),
						color = color))

		for zi in range(num_z_levels):
			z = z_levels[zi]
			r_inner = radius*( z - 0.5*dz)
			r_outer = radius*( z + 0.5*dz)
			theta_angle = 360* theta_levels[ti]/(2*np.pi )

			##### MY map!
			if with_deuteranomaly == False:
				color = sRGB1_radial_image[ti,zi]
			else:
				color =  deut_radial_image[ti,zi]

			ax.add_patch(mpatches.Wedge(center = center, r=
			r_outer, width= r_outer-r_inner,
				theta1 = theta_angle-360*dtheta/2/(2*np.pi ), theta2 = theta_angle+360*dtheta/2/(2*np.pi ) ,
				color = color))

	label = 'Colorwheel'
	if with_deuteranomaly:
		label = label + ' with\nDeuteranomaly'

	ax.text(xc, yc + radius*1.1, label , va='bottom', ha='center')
	ax.axis('equal')

def plot_small_wave_angle_test(ax, my_map, nwaves = 30, amplitude = 0.2, pixels_per_wave = 8):

	n_theta_waves = nwaves

	num_z_levels =  64
	num_theta_levels = pixels_per_wave*(n_theta_waves+3) # a few extra pixels help

	dtheta = 1.0/num_theta_levels
	theta_levels = np.arange(dtheta/2.0, 1.0, dtheta)*np.pi*2
	dtheta = dtheta*np.pi*2
	dz = 1.0/num_z_levels
	z_levels = np.arange(dz/2.0, 1.0, dz)

	line_len = 0.08
	### hue
	theta_grid_2D, z_grid_2D =  np.meshgrid(theta_levels, z_levels, sparse=False, indexing='xy')
	theta_grid_2D = theta_grid_2D + z_grid_2D*amplitude*np.sin(n_theta_waves*theta_grid_2D)

	cvd_space = {"name": "sRGB1+CVD",
		"cvd_type": "deuteranomaly",
		"severity": 50}


	image = my_map(theta_grid_2D, np.ones(z_grid_2D.shape), verbose = False)
	deut_image = np.clip(cspace_convert(image, cvd_space, "sRGB1"),0,1)

	ax.imshow(image,      origin = 'lower', extent = (0, 2*np.pi, 1,  2), interpolation ='bilinear' )
	ax.imshow(deut_image, origin = 'lower', extent = (0, 2*np.pi, 0,  1), interpolation ='bilinear' )
	ax.vlines(my_map.theta_knots, [0], [line_len], color = 'k', linewidth = 0.75)
	ax.vlines(my_map.theta_knots, [1], [1+line_len], color = 'k', linewidth = 0.75)
	ax.set_ylim(0,2*np.pi)
	ax.axis('auto')

def plot_small_wave_magnitude_test(ax, my_map, nwaves = 30, amplitude = 0.1, pixels_per_wave = 8):

	n_z_waves = nwaves

	num_theta_levels = 128
	num_z_levels =     pixels_per_wave*(n_z_waves+3)

	dtheta = 1.0/num_theta_levels
	theta_levels = np.arange(dtheta/2.0, 1.0, dtheta)*np.pi*2
	dtheta = dtheta*np.pi*2
	dz = 1.0/num_z_levels
	z_levels = np.arange(dz/2.0, 1.0, dz)

	line_len = 0.08
	### brightness

	theta_grid_2D, z_grid_2D =  np.meshgrid(theta_levels, z_levels, sparse=False, indexing='xy')
	z_grid_2D =      z_grid_2D     + theta_grid_2D/(2*np.pi)*amplitude*np.sin(n_z_waves*2*np.pi*z_grid_2D)

	cvd_space = {"name": "sRGB1+CVD",
		"cvd_type": "deuteranomaly",
		"severity": 50}


	image = my_map(theta_grid_2D, z_grid_2D,verbose = False)
	deut_image = np.clip(cspace_convert(image, cvd_space, "sRGB1"),0,1)

	ax.imshow(image,origin = 'lower',      extent = (0, 1, 0.0, 2*np.pi) , interpolation ='bilinear' )
	ax.imshow(deut_image,origin = 'lower', extent = (1, 2, 0.0, 2*np.pi), interpolation ='bilinear' )

	ax.hlines(my_map.theta_knots,[1-line_len/2.0], [1+line_len/2.0] ,   color = 'w', linewidth = 1.5)
	ax.hlines(my_map.theta_knots,[1-line_len/2.0], [1+line_len/2.0],   color = 'k', linewidth = 0.8)
	#ax.hlines([1-line_len/2.0], [1], my_map.theta_knots+1, color = 'k', linewidth = 0.75)

	ax.set_ylim(0,2*np.pi)
	ax.set_xlim(0,2)
	ax.axis('auto')



def plot_analytic_vortex_test_data(ax, my_map, npts = 16, xmax = 1.5, ymax = 1.75, with_quiver = True, with_deuteranomaly = False):

	X, Y = np.meshgrid(linspace(0, xmax , xmax * npts), linspace(0, ymax , ymax*npts))
	U = np.sin(2*pi*Y)
	V = np.cos(2*pi*X)

	angle = arctan2(V,U)
	magnitude = sqrt(V**2 + U**2)

	magnitude_norm = magnitude/magnitude.max()
	angle_norm     = angle

	image = my_map(angle_norm, magnitude_norm)
	if with_deuteranomaly:
		cvd_space = {"name": "sRGB1+CVD",
				"cvd_type": "deuteranomaly",
				"severity": 50}
		image =clip( cspace_convert(image, cvd_space, "sRGB1"),  0,1)

	ax.imshow(image, origin = 'lower',      extent = (0, xmax, 0, ymax) , interpolation ='None' )
	if with_quiver: ax.quiver(X, Y, U, V, units='width')



def cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False, color_space = default_colorspace, ab_grid = 512):



	print('\n-->testing %s\n'%name)

	my_map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = color_space)
	#my_map.min_theta_pts

	my_map.name = name

	from matplotlib.pylab import subplots, show, subplot
	nrows = 3
	ncols = 3
	subsize = 3.0
	fig, axes =  subplots(nrows=nrows,  ncols=ncols, figsize=(subsize*ncols, subsize*nrows))
	fig.suptitle(my_map.name, fontsize = 14)


	plot_knots_on_isoluminant_slice(axes[0,0], my_map, ab_grid)

	plot_arc_length_vs_spline_parameter_t(axes[0,1], my_map)
	plot_perceptual_derivative(axes[1,0], my_map,)
	plot_perceptual_derivative(axes[1,1], my_map, use_delta_instead = True)


	plot_colorwheel(axes[2,1], my_map)
	plot_colorwheel(axes[2,1], my_map,
				num_theta_levels = 32,
				num_z_levels = 16,
				center = (3,0),
				radius=1.0,
				with_side_bar = False,
				with_deuteranomaly = True)

	plot_small_wave_angle_test(axes[0,2], my_map)
	plot_small_wave_magnitude_test(axes[1,2], my_map)

	if False:
		plot_analytic_vortex_test_data(axes[2,0], my_map)
		axes[2,0].set_title('This Map')

		plot_analytic_vortex_test_data(axes[2,2], my_map, with_deuteranomaly = True)
		axes[2,2].set_title('Deuteranomaly')

	else:
		########### sample data!
		z_array = np.loadtxt('example_data/magnitudes.txt')
		theta_array = np.loadtxt('example_data/theta_angles.txt')

		image = my_map(theta_array*2*np.pi, z_array)
		cvd_space = {"name": "sRGB1+CVD",
				"cvd_type": "deuteranomaly",
				"severity": 50}
		deut_image = np.clip( cspace_convert(image, cvd_space, "sRGB1"),  0,1)


		axes[2,0].imshow(image, origin = 'lower')
		axes[2,0].set_title('This Map')

		axes[2,2].imshow(deut_image, origin = 'lower')
		axes[2,2].set_title('Deuteranomaly')


	fig.tight_layout(pad =0.1)
	fig.subplots_adjust(top = 0.95)
	fig.savefig(my_map.name+'_colormap.pdf',transparent = True, dpi =1200)
	fig.savefig(my_map.name+'_colormap.svg',transparent = True, dpi =1200)
	########
	if show_output:show()


if __name__ == "__main__":
	from colorspacious import cspace_convert
	from colorsys import hsv_to_rgb

	lab_color = cspace_convert(hsv_to_rgb(0.0,1.0,1.0),"sRGB1",  "CAM02-UCS")
	HSV_ab_angle0 = np.arctan2(lab_color[2], lab_color[1])
	print(lab_color, HSV_ab_angle0 * 180/pi)


	## for CAM02-UCS
	L_max = 73
	radius = 26
	name = 'circle-%.1f-%.1f' %(radius,L_max)
	center = (0,0)
	theta = np.arange(0, 2*pi, 2*pi/8)
	a_knots = center[0] + radius*cos(theta + HSV_ab_angle0)
	b_knots = center[1] + radius*sin(theta + HSV_ab_angle0)

	cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = True)

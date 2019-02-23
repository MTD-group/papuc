

#######################
from numpy import array, pi, linspace, mod,  append, arange, sqrt, diff, mean, clip, zeros, any, all, cos, sin, floor, meshgrid, ones
from scipy.interpolate import splrep, splev



from tools import periodic_spline, periodic_spline_test, isoluminant_uniform_spline_colormap
from tools import HSV_colormap

if False:
	periodic_spline_test()







def cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = False, lab_grid = 1024):

	na_grid = nb_grid = lab_grid
	num_theta_levels = 32
	num_z_levels = 16
	print('\n-->testing %s\n'%name)

	my_map = isoluminant_uniform_spline_colormap(a_knots, b_knots, L_max, color_space = 'CIELab')
	my_map.min_theta_pts



	####### HSV comparison
	from colorspacious import cspace_convert
	from colorsys import hsv_to_rgb
	HSV_pts = 300
	HSV_V_levels = [0.2, 0.5, 1.0]
	L_HSV = zeros((len(HSV_V_levels),HSV_pts))
	a_HSV = zeros((len(HSV_V_levels),HSV_pts))
	b_HSV = zeros((len(HSV_V_levels),HSV_pts))
	sRGB1_HSV = zeros((len(HSV_V_levels),HSV_pts,3))
	v_knots_HSV_fine, dv_HSV_fine = linspace(0,1,HSV_pts, retstep=True)
	v_knots_HSV_fine_centers = linspace(0,1,HSV_pts-1)

	for HSV_V_index in range(len(HSV_V_levels)):
		HSV_V_level = HSV_V_levels[HSV_V_index]
		for vi in range(len(v_knots_HSV_fine)):
			v = v_knots_HSV_fine[vi]
			sRGB1_HSV[HSV_V_index,vi] = hsv_to_rgb( v , 1.0, HSV_V_level)
			lab_color = cspace_convert(sRGB1_HSV[HSV_V_index, vi], "sRGB1", my_map.color_space )
			L_HSV[HSV_V_index,vi]= lab_color[0]
			a_HSV[HSV_V_index,vi]= lab_color[1]
			b_HSV[HSV_V_index,vi]= lab_color[2]

	delta_L = diff(L_HSV[-1])
	delta_a = diff(a_HSV[-1])
	delta_b = diff(b_HSV[-1])

	dE_HSV_fine = sqrt(delta_L**2 + delta_a**2 + delta_b**2)
	#dE_HSV_fine = sqrt( delta_a**2 + delta_b**2)


	######################



	import matplotlib.patches as mpatches
	from matplotlib.collections import PatchCollection


	color_grid = True

	from matplotlib.pylab import subplots, show, subplot
	nrows = 2
	ncols = 2
	subsize = 3.0
	fig1, axes1 =  subplots(nrows=1,      ncols=ncols, figsize=(subsize*ncols, subsize*1))
	fig2, axes2 =  subplots(nrows=nrows , ncols=ncols, figsize=(subsize*ncols, subsize*nrows))
	fig_small_wave, axes_small_wave =  subplots(nrows=nrows , ncols=ncols, figsize=(subsize*ncols, subsize*nrows))
	fig_img, axes_img =  subplots(nrows=nrows , ncols=ncols, figsize=(subsize*ncols, subsize*nrows))

	axes1[0].set_title("Color Map '%s'"%name)

	axes1[0].plot(a_HSV[-1], b_HSV[-1], color = 'lightgrey')
	### plot the smooth spline path
	theta_plot = linspace(0,1.0,my_map.min_theta_pts*4)
	axes1[0].plot(
			my_map.a_func_of_theta(theta_plot),
			my_map.b_func_of_theta(theta_plot), color ='k')

	### plot the equally spaced points from the minimum v_refinement
	theta_max_space = linspace(0,1.0,my_map.min_theta_pts)
	axes1[0].plot(
			my_map.a_func_of_theta(theta_max_space),
			my_map.b_func_of_theta(theta_max_space), marker = 'x', color = 'k',linestyle = '' )
	### plot the spline knots
	axes1[0].plot(my_map.a_knots,my_map.b_knots, marker = 'o', linestyle = '', mfc = 'white', mec = 'black', markersize = 10)
	for i in range(len(my_map.a_knots)):
		axes1[0].text(
					my_map.a_knots[i],
					my_map.b_knots[i], '%i'%i, ha = 'center', va = 'center', fontsize = 8)

	axes1[0].set_xlabel('a*')
	axes1[0].set_ylabel('b*')
	axes1[0].axis('equal')

	if color_grid:
		#na_grid = nb_grid = 256# 150

		#amin, amax = a_fine_for_plot.min()*1.1, a_fine_for_plot.max()*1.1
		#bmin, bmax = b_fine_for_plot.min()*1.1, b_fine_for_plot.max()*1.1

		# good for L*=75
		#amin, amax = -80, 65
		#bmin, bmax = -45, 85

		amin, amax = -100, 100
		bmin, bmax = -100, 100





		use_image = True
		if use_image:
			print('making image')
			da = (amax-amin)/na_grid
			db = (bmax-bmin)/nb_grid
			a_points = linspace(amin+da/2, amax-da/2, na_grid)
			b_points = linspace(bmin+db/2, bmax-db/2, nb_grid)

			a_grid, b_grid = meshgrid(a_points,b_points, indexing = 'xy')

			lab_map = ones((na_grid, nb_grid, 3))*L_max

			lab_map[:,:,1] = a_grid
			lab_map[:,:,2] = b_grid

			bg_color = array([1, 1, 1, 0])

			sRGB1_map = zeros((na_grid, nb_grid, 4)) # default is transparent
			sRGB1_map[:,:,0:3] = cspace_convert(lab_map, my_map.color_space , "sRGB1")
			for a_index in range(na_grid):
				for b_index in range(nb_grid):
					sRGB1_color = sRGB1_map[a_index,b_index,0:3]
					if  all(0<sRGB1_color) and all(sRGB1_color<1):
						sRGB1_map[a_index,b_index,3] = 1
					else:
						sRGB1_map[a_index,b_index] = bg_color

			axes1[0].imshow(sRGB1_map, extent=(amin,amax,bmin,bmax),
				 origin = 'lower', interpolation = 'bicubic') #'none', 'nearest', 'bilinear', 'bicubic'

		else: # This will use slow patches
			a_edges, da = linspace(amin,amax, na_grid+1, retstep=True)
			b_edges, db = linspace(bmin,bmax, nb_grid+1, retstep=True)

			sRGB1_edgecolor = cspace_convert([100-L_max,0,0], my_map.color_space , "sRGB1")

			for a_index in range(na_grid):
				for b_index in range(nb_grid):
					a = a_edges[a_index]+0.5*da
					b = b_edges[b_index]+0.5*db

					lab_color = [L_max, a, b]
					sRGB1_color = cspace_convert(lab_color, my_map.color_space , "sRGB1")

					if  all(0<sRGB1_color) and all(sRGB1_color<1):

						axes1[0].add_patch(mpatches.Rectangle(
								xy = (a_edges[a_index], b_edges[b_index]),
								width = da,
								height = db,
								color = sRGB1_color))#,
								#linestyle = '-',linewidth = 4.0,
								#edgecolor = 'k'))


	#############
	####### arc length
	arc_len_in_t =linspace(0.0,1.0, my_map.min_t_pts)
	arc_len = my_map.arc_length * my_map.func_theta_of_t(arc_len_in_t)
	arc_len[-1] = my_map.arc_length
	axes2[0,0].plot(arc_len_in_t, arc_len )

	t_at_theta_max_space = my_map.func_t_of_theta(theta_max_space)
	t_at_theta_max_space[-1] = 1.0
	axes2[0,0].plot(t_at_theta_max_space, theta_max_space* my_map.arc_length, marker = '+', linestyle = '', color = 'k')
	axes2[0,0].set_xlabel('Spline Parametric, t')
	axes2[0,0].set_ylabel('Arc Length')
	axes2[0,0].set_xlim((0,1))
	################
	# these are pairwise differences
	#t_knots = arange(0,1.0, 1.0/len(a_knots))
	#v_knots = func_v_of_t(t_knots)
	for theta in my_map.theta_knots:
		axes2[0,1].axvline(theta, color = 'grey')

	dE = my_map.compute_min_theta_points_for_uniformity()
	dv = 1.0/(len(dE))
	v_points_for_center_diff = arange(dv/2.0, 1.0, dv)
	axes2[0,1].plot( v_points_for_center_diff, dE, marker = 'o', markersize = 2, label = '%i pts'%(len(dE)+1))
	axes2[0,1].set_xlim((0,1))
	axes2[0,1].legend(handlelength = 1.0)
	axes2[0,1].set_xlabel(r'Mapped Phase $\theta$,'+'\nNormalized to Arc length')
	axes2[0,1].set_ylabel('$\Delta E$')
	###############

	#axes[0,3] = subplot(nrows = 3, ncols = 4, index = 3, projection='polar')
	#subplot(axes[0,3], projection='polar')


	###################
	axes2[1,0].plot( v_points_for_center_diff, dE/dv, marker = 'o', markersize = 2, label = 'This Map, z = 1.0')#label = '%i pts'%(len(dE)+1))
	axes2[1,0].plot(v_knots_HSV_fine_centers, dE_HSV_fine/dv_HSV_fine, markersize = 2, label = r'HSV, V = %.1f'% HSV_V_level)
	axes2[1,0].set_xlim((0,1))
	axes2[1,0].legend(handlelength = 1.0)
	axes2[1,0].set_xlabel(r'Mapped Phase $\theta$,'+'\nNormalized to Arc length')
	axes2[1,0].set_ylabel(r'$\frac{dE}{dv}$')


	#################

	axes2[1,1].axhline(L_max, label = 'This Map, z = 1.0')#label = '%i pts'%(len(dE)+1))


	for HSV_V_index in [-1]:#range(len(HSV_V_levels)-1,-1,-1):
			HSV_V_level = HSV_V_levels[HSV_V_index]
			axes2[1,1].plot(v_knots_HSV_fine, L_HSV[HSV_V_index], color = 'C1', label = r'HSV, V = %.1f'% HSV_V_level)
	axes2[1,1].set_xlim((0,1))
	axes2[1,1].set_ylim((0,100))
	axes2[1,1].legend(handlelength = 1.0)
	axes2[1,1].set_xlabel(r'Mapped Phase, $\theta$')
	axes2[1,1].set_ylabel(r'$L^{*}$')



	############ make a color map


	############
	dtheta = 1.0/num_theta_levels
	theta_levels = arange(dtheta/2.0, 1.0, dtheta)
	dz = 1.0/num_z_levels
	z_levels = arange(dz/2.0, 1.0, dz)
	#################

	cvd_space = {"name": "sRGB1+CVD",
		"cvd_type": "deuteranomaly",
		"severity": 50}


	theta_grid_2D, z_grid_2D =  meshgrid(theta_levels, z_levels, sparse=False, indexing='ij')

	sequence = my_map(theta_levels, ones(theta_levels.shape))
	deut_sequence = cspace_convert(sequence, cvd_space, "sRGB1")
	image = my_map(theta_grid_2D, z_grid_2D, verbose = False)
	deut_image = cspace_convert(image, cvd_space, "sRGB1")

	HSV_sequence = HSV_colormap(theta_levels, ones(theta_levels.shape))
	HSV_deut_sequence = cspace_convert(HSV_sequence, cvd_space, "sRGB1")
	HSV_image = HSV_colormap(theta_grid_2D, z_grid_2D, verbose = False)
	HSV_deut_image = cspace_convert(HSV_image, cvd_space, "sRGB1")

	my_map_y = 1.1
	HSV_y = -my_map_y

	newax=axes1[1]
	for i in range(num_theta_levels):

		### MY map!
		sRGB1_color = sequence[i]
		deut_color =  deut_sequence[i]
		newax.add_patch(mpatches.Rectangle(
					xy = (-0.4,  2.0*(theta_levels[i]-dtheta/2.0)-1.0 + my_map_y),
					width = 0.2,
					height = dtheta*2,
					color = sRGB1_color))
		# deuteranomaly version
		deut_color = clip(deut_color, 0,1)
		newax.add_patch(mpatches.Rectangle(
					xy = (2.2,  2.0*(theta_levels[i]-dtheta/2.0)-1.0+ my_map_y),
					width = 0.2,
					height = dtheta*2,
					color = deut_color))

		#####HSV version
		sRGB1_color = HSV_sequence[i]
		deut_color =  HSV_deut_sequence[i]
		newax.add_patch(mpatches.Rectangle(
					xy = (-0.4,  2.0*(theta_levels[i]-dtheta/2.0)-1.0 + HSV_y),
					width = 0.2,
					height = dtheta*2,
					color = sRGB1_color))
		# deuteranomaly version
		deut_color = clip(deut_color, 0,1)
		newax.add_patch(mpatches.Rectangle(
					xy = (2.2,  2.0*(theta_levels[i]-dtheta/2.0)-1.0 + HSV_y),
					width = 0.2,
					height = dtheta*2,
					color = deut_color))

		for zi in range(num_z_levels):
			z = z_levels[zi]
			r_inner = z-dz/2.0
			r_outer = z+dz/2.0
			theta_angle = 360* theta_levels[i]

			##### MY map!
			sRGB1_color = image[i,zi]
			deut_color =  deut_image[i,zi]
			newax.add_patch(mpatches.Wedge(center = (-1.5, my_map_y), r=
			r_outer, width= r_outer-r_inner,
				theta1 = theta_angle-360*dtheta/2, theta2 = theta_angle+360*dtheta/2 ,
				color = sRGB1_color))
			# deuteranomaly version
			deut_color = clip(deut_color, 0,1)
			newax.add_patch(mpatches.Wedge(center = (1,  my_map_y), r=
			r_outer, width= r_outer-r_inner,
				theta1 = theta_angle-360*dtheta/2, theta2 = theta_angle+360*dtheta/2 ,
				color = deut_color))

			##### HSV
			sRGB1_color = HSV_image[i,zi]
			deut_color =  HSV_deut_image[i,zi]
			newax.add_patch(mpatches.Wedge(center = (-1.5, HSV_y), r=
			r_outer, width= r_outer-r_inner,
				theta1 = theta_angle-360*dtheta/2, theta2 = theta_angle+360*dtheta/2 ,
				color = sRGB1_color))
			# deuteranomaly version
			deut_color = clip(deut_color, 0,1)
			newax.add_patch(mpatches.Wedge(center = (1, HSV_y), r=
			r_outer, width= r_outer-r_inner,
				theta1 = theta_angle-360*dtheta/2, theta2 = theta_angle+360*dtheta/2 ,
				color = deut_color))

	newax.text(-1.5, 1.1+my_map_y, 'Colorwheel', va='bottom', ha='center')
	newax.text(1.2, 1.1+my_map_y, 'Deuteranomaly', va='bottom', ha='center')

	newax.set_ylim(-1,1)
	newax.set_xlim(-1,1)
	newax.axis('equal')
	### HSV
	newax.text(-1.5, -1.1+HSV_y, 'HSV\nColorwheel', va='top', ha='center')
	newax.text(1.2, -1.1+HSV_y, 'HSV\nDeuteranomaly', va='top', ha='center')

	################# small wave testing of map


	n_z_waves = n_theta_waves = 40

	amplitude = 0.1

	num_theta_levels = 5*n_theta_waves+20
	num_z_levels =   5*n_z_waves+20
	dtheta = 1.0/num_theta_levels
	theta_levels = arange(dtheta/2.0, 1.0, dtheta)
	dz = 1.0/num_z_levels
	z_levels = arange(dz/2.0, 1.0, dz)

	line_len = 0.08
	### hue
	theta_grid_2D, z_grid_2D =  meshgrid(theta_levels, z_levels, sparse=False, indexing='xy')
	theta_grid_2D = theta_grid_2D + z_grid_2D*amplitude*sin(n_theta_waves*2*pi*theta_grid_2D)

	image = my_map(theta_grid_2D, ones(z_grid_2D.shape), verbose = False)
	deut_image = clip(cspace_convert(image, cvd_space, "sRGB1"),0,1)

	axes_small_wave[0,0].imshow(image,origin = 'lower', extent =      (0,  1, 1,  2), interpolation ='bilinear' )
	axes_small_wave[0,0].imshow(deut_image,origin = 'lower', extent = (0, 1, 0,   1), interpolation ='bilinear' )
	axes_small_wave[0,0].set_ylim(0,1)
	axes_small_wave[0,0].axis('auto')

	axes_small_wave[0,0].vlines(my_map.theta_knots, [0], [line_len], color = 'k', linewidth = 0.75)
	axes_small_wave[0,0].vlines(my_map.theta_knots, [1], [1+line_len], color = 'k', linewidth = 0.75)

	########brightness


	theta_grid_2D, z_grid_2D =  meshgrid(theta_levels, z_levels, sparse=False, indexing='xy')
	z_grid_2D =      z_grid_2D     + theta_grid_2D*amplitude*sin(n_z_waves*2*pi*z_grid_2D)

	image = my_map(theta_grid_2D, z_grid_2D,verbose = False)
	deut_image = clip(cspace_convert(image, cvd_space, "sRGB1"),0,1)

	axes_small_wave[0,1].imshow(image,origin = 'lower',      extent = (0, 1, 0.0, 1.0) , interpolation ='bilinear' )
	axes_small_wave[0,1].imshow(deut_image,origin = 'lower', extent = (1, 2, 0.0, 1.0), interpolation ='bilinear' )
	axes_small_wave[0,1].set_xlim(0,1)
	axes_small_wave[0,1].axis('auto')


	axes_small_wave[0,1].vlines(my_map.theta_knots, [1-line_len/2.0], [1], color = 'k', linewidth = 0.75)
	axes_small_wave[0,1].vlines(my_map.theta_knots+1, [1-line_len/2.0], [1], color = 'k', linewidth = 0.75)


	# ##############small wave testing of HSV

	### hue
	theta_grid_2D, z_grid_2D =  meshgrid(theta_levels, z_levels, sparse=False, indexing='xy')
	theta_grid_2D = theta_grid_2D + z_grid_2D*amplitude*sin(n_theta_waves*2*pi*theta_grid_2D)

	image = HSV_colormap(theta_grid_2D, ones(z_grid_2D.shape), verbose = False)
	deut_image = clip(cspace_convert(image, cvd_space, "sRGB1"),0,1)

	axes_small_wave[1,1].imshow(image,origin = 'lower', extent =      (0,  1, 1,  2), interpolation ='bilinear' )
	axes_small_wave[1,1].imshow(deut_image,origin = 'lower', extent = (0, 1, 0,   1), interpolation ='bilinear' )
	axes_small_wave[1,1].set_ylim(0,2)
	axes_small_wave[1,1].axis('auto')

	########brightness


	theta_grid_2D, z_grid_2D =  meshgrid(theta_levels, z_levels, sparse=False, indexing='xy')
	z_grid_2D =      z_grid_2D     + theta_grid_2D*amplitude*sin(n_z_waves*2*pi*z_grid_2D)

	image = HSV_colormap(theta_grid_2D, z_grid_2D,verbose = False)
	deut_image = clip(cspace_convert(image, cvd_space, "sRGB1"),0,1)

	axes_small_wave[1,0].imshow(image,origin = 'lower',      extent = (0.0, 1, 0.0, 1.0) , interpolation ='bilinear' )
	axes_small_wave[1,0].imshow(deut_image,origin = 'lower', extent = (1, 2, 0.0, 1.0), interpolation ='bilinear' )
	axes_small_wave[1,0].set_xlim(0,2)
	axes_small_wave[1,0].axis('auto')




	########### sample data!
	from numpy import loadtxt

	z_array = loadtxt('example_data/magnitudes.txt')
	theta_array = loadtxt('example_data/theta_angles.txt')

	image = my_map(theta_array, z_array)
	deut_image =clip( cspace_convert(image, cvd_space, "sRGB1"),  0,1)


	axes_img[0,0].imshow(image, origin = 'lower')
	axes_img[0,0].set_title('This Map')

	cvd_space = {"name": "sRGB1+CVD",
		"cvd_type": "deuteranomaly",
		"severity": 50}
	axes_img[1,0].imshow(deut_image, origin = 'lower')
	axes_img[1,0].set_title('Deuteranomaly')



	HSV_image = HSV_colormap(mod(theta_array,1.0), z_array)
	axes_img[0,1].imshow(HSV_image, origin = 'lower')
	axes_img[0,1].set_title('HSV')

	deut_HSV_image =clip( cspace_convert(HSV_image, cvd_space, "sRGB1"),  0,1)
	axes_img[1,1].imshow(deut_HSV_image, origin = 'lower')
	axes_img[1,1].set_title('HSV Deuteranomaly')

	#######
	#print('Saving pdf')
	#fig.tight_layout(pad =0.1)
	#fig.subplots_adjust(wspace = 0.32)
	#fig.savefig(name+'_colormap.pdf',transparent = True, dpi =1200)
	#fig.savefig(name+'_colormap.svg',transparent = True, dpi =1200)
	########
	if show_output:show()




if __name__ == "__main__":
	from colorspacious import cspace_convert
	from colorsys import hsv_to_rgb
	from numpy import arctan2
	lab_color = cspace_convert(hsv_to_rgb(0.0,1.0,1.0),"sRGB1",  "CIELab")
	HSV_ab_angle0 = arctan2(lab_color[2], lab_color[1])
	print(lab_color, HSV_ab_angle0 * 180/pi)


	## for CIELab
	L_max = 74
	radius = 40
	name = 'circle-%.1f-%.1f' %(radius,L_max)
	center = (0,0)
	theta = arange(0, 2*pi, 2*pi/8)
	a_knots = center[0] + radius*cos(theta + HSV_ab_angle0)
	b_knots = center[1] + radius*sin(theta + HSV_ab_angle0)

	cyclic_colormap_test(a_knots, b_knots, L_max, name, show_output = True, lab_grid =512)

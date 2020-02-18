


import numpy as np
from scipy.interpolate import splrep, splev
from colorspacious import cspace_convert


default_colorspace = "CAM02-UCS"

def periodic_spline(theta_knots, f_knots, period = 2*np.pi):
	### making some local copies ###
	theta_knots_local = np.array(theta_knots)
	f_knots_local     = np.array(f_knots)

	theta0 = theta_knots_local[0]

	### takes care of non_explicit periodicity oddness
	theta_knots_local =  np.append( theta_knots_local, [theta0 + period] )
	f_knots_local =  np.append (f_knots_local, [f_knots_local[0]+period ]  )

	tck = splrep(theta_knots_local, f_knots_local,
			xb = theta0, xe = theta0 + period,
			 per = True, s=0 )

	def my_func(theta):
		#values = splev(theta, tck) # gotta have the modulo
		values = splev( np.mod(theta-theta0 ,period) +theta0, tck)
		return values

	def my_func_prime(theta):
		#values = splev(theta, tck) # gotta have the modulo
		values = splev( np.mod(theta-theta0 ,period) +theta0, tck, der= 1)
		return values

	return my_func, my_func_prime



def periodic_spline_test():
	a_knots =  [ -30, 30, 40, -15,   -60]
	b_knots =  [  -20, -35 , -10, 15,  40]

	theta_knots  = np.array([ 45, 125, 230, 320, 340 ])

	a_func, a_func_prime = periodic_spline(theta_knots, a_knots, period = 360)
	b_func, b_func_prime = periodic_spline(theta_knots, b_knots, period = 360)

	#print(a_func(theta_knots))

	from matplotlib.pylab import plot, axvline, gca, show
	from matplotlib.ticker import MultipleLocator

	theta_smooth =  np.arange(-45, 360+45,1)

	plot(theta_smooth, a_func(theta_smooth), color = 'teal')
	plot(theta_smooth, b_func(theta_smooth), color = 'orange')

	plot(theta_knots,a_knots, marker = 'o', linestyle = '', color = 'teal')
	plot(theta_knots,b_knots, marker = 'o', linestyle = '', color = 'orange')

	axvline([0],color = 'k',linestyle = '--')
	axvline([360],color = 'k',linestyle = '--')

	gca().set_xlim(-45,360+45)
	gca().xaxis.set_major_locator(MultipleLocator(45))
	gca().set_xlabel('cycle degrees')
	show()




class isoluminant_uniform_spline_colormap:
	def __init__(self,a_knots, b_knots, L_max,
		delta_E_model = '',
		abs_cumulative_error_cutoff = 0.01, initial_arclen_npts = 4, color_space = default_colorspace ):

		self.color_space = color_space # you can also use "CIELab"

		#abs_cumulative_error_cutoff = 0.01 # about 1% of what a human could even detect!
		#initial_arclen_npts = 4 pretty coarse

		#print('\n-->testing %s\n'%name)
		t_knots =  np.arange(0,1.0, 1.0/len(a_knots))
		a_func, a_func_prime = periodic_spline(t_knots, a_knots, period = 1.0 )
		b_func, b_func_prime = periodic_spline(t_knots, b_knots, period = 1.0 )

		def differtial_arc(t):
			return  np.sqrt(a_func_prime(t)**2 + b_func_prime(t)**2)
		from scipy.integrate import simps, cumtrapz
		#print(theta_knots)

		# you only cared about arc length this a good metric, we need the integrand
		#arc_len_quad, quad_error = quad(differtial_arc, 0, 1)

		t_pts = initial_arclen_npts
		t_knots_fine, dt = np.linspace(0, 1.0, t_pts, retstep=True)
		cum_int  = cumtrapz(differtial_arc(t_knots_fine), dx = dt, initial = 0.0)
		#abs_error = abs(cum_int[-1]- arc_len_quad)
		abs_error = abs_cumulative_error_cutoff*2


		while abs_error > abs_cumulative_error_cutoff:
			cum_int_old = cum_int
			t_pts = t_pts*2
			t_knots_fine, dt = np.linspace(0, 1.0, t_pts, retstep=True)
			cum_int  = cumtrapz(differtial_arc(t_knots_fine),dx = dt, initial = 0.0)
			abs_error = abs(cum_int[-1] - cum_int_old[-1])

			print('Pathlength Error %0.3e t_pts %i'%( abs_error, t_pts))

		self.min_t_pts = t_pts

		print('Previous Arc Length %0.3e Current Arc Length %0.3e t_pts %i'%( cum_int_old[-1], cum_int[-1],  t_pts))

		self.arc_length = cum_int[-1]

		print('arc len',self.arc_length)
		theta_knots_fine =  cum_int/(self.arc_length)*2*np.pi

		from scipy.interpolate import interp1d

		func_theta_of_t_raw = interp1d(
			x = t_knots_fine,
			y = theta_knots_fine,
			kind = 'cubic')
		def func_theta_of_t(t):
			# this fixes the looping issue
			return func_theta_of_t_raw(t-np.floor(t))
		self.func_theta_of_t = func_theta_of_t


		func_t_of_theta_raw = interp1d(
			x = theta_knots_fine,
			y = t_knots_fine,
			kind = 'cubic')
		def func_t_of_theta(theta):
			# this fixes the looping issue
			return func_t_of_theta_raw( theta-2*np.pi*np.floor(theta/(2*np.pi)))
		self.func_t_of_theta = func_t_of_theta



		def a_func_of_theta(theta):
			return  a_func(func_t_of_theta(theta))
		self.a_func_of_theta = a_func_of_theta

		def b_func_of_theta(theta):
			return  b_func(func_t_of_theta(theta))
		self.b_func_of_theta = b_func_of_theta

		self.L_max = L_max
		self.a_knots = a_knots
		self.b_knots = b_knots
		self.theta_knots = func_theta_of_t(t_knots)

		self.compute_min_theta_points_for_uniformity()


	def compute_min_theta_points_for_uniformity(self, JND = None, max_delta_E_deviation = 1.0):

		if JND == None:
			max_delta_E = 10.0
		else:
			max_delta_E = JND
		###### how many theta_knots do we need? this computes it
		theta_pts =  len(self.a_knots)*2+1
		theta_knots_fine, dtheta = np.linspace(0, 1.0, theta_pts, retstep=True)
		delta_a =  np.diff(self.a_func_of_theta(theta_knots_fine))
		delta_b =  np.diff(self.b_func_of_theta(theta_knots_fine))
		dE =  np.sqrt(delta_a**2 + delta_b**2)# scalling by 1/dtheta was ugly
		dE_list = [dE]
		mean_dE =    np.mean(dE)
		max_above = dE.max()-mean_dE
		min_below = mean_dE - dE.min()
		#print('mean_delta_E', mean_dE, 'max_above', max_above, 'min_below', min_below, 'uniform_maping theta_pts', theta_pts  )
		print('mean_delta_E %0.3e max_above %0.3e min_below %0.3e uniform_maping theta_pts %i'%( mean_dE, max_above, min_below, theta_pts  ))


		while max(max_above,min_below) > max_delta_E_deviation or dE.max() > max_delta_E:
			theta_pts = theta_pts*2
			theta_knots_fine, dtheta = np.linspace(0, 1.0, theta_pts, retstep=True)
			delta_a =  np.diff(self.a_func_of_theta(theta_knots_fine))
			delta_b =  np.diff(self.b_func_of_theta(theta_knots_fine))
			dE =  np.sqrt(delta_a**2 + delta_b**2)
			dE_list.append(dE)
			mean_dE =    np.mean(dE)
			max_above = dE.max()-mean_dE
			min_below = mean_dE - dE.min()
			print('mean_delta_E %0.3e max_above %0.3e min_below %0.3e uniform_maping theta_pts %i'%( mean_dE, max_above, min_below, theta_pts  ))
			#print('mean_delta_E', mean_dE, 'max_above', max_above, 'min_below', min_below, 'uniform_maping theta_pts', theta_pts  )

		self.min_theta_pts = theta_pts
		return dE

	def __call__(self, theta, z, scale_alpha = False, clip_values = True, verbose = True):


		z_ar = np.array(z)
		theta_ar =  np.array(theta)

		dims = len(z_ar.shape)

		if scale_alpha:
			image = np.zeros(list(z_ar.shape)+[4])
		else:
			image = np.zeros(list(z_ar.shape)+[3])

		if dims == 1:
			for i in range(z_ar.shape[0]):
				a = self.a_func_of_theta(theta_ar[i])
				b = self.b_func_of_theta(theta_ar[i])

				lab_color = z_ar[i]*np.array([self.L_max,a,b])
				if scale_alpha: lab_color = np.array([self.L_max,a,b])

				sRGB1_color = cspace_convert(lab_color, self.color_space, "sRGB1")
				if clip_values:
					if any(sRGB1_color>1) or any(sRGB1_color<0):
						sRGB1_color_clipped =  np.clip( sRGB1_color, 0,1)
						Lab_color_clipped = cspace_convert(sRGB1_color_clipped, "sRGB1", self.color_space)
						if verbose:print (self.color_space+'/sRGB1 color', np.array(lab_color), '/', sRGB1_color,
							'\n-----> clipped to',  Lab_color_clipped, '/', sRGB1_color_clipped)
						sRGB1_color = sRGB1_color_clipped

				if scale_alpha: sRGB1_color =  np.append(sRGB1_color, z_ar[i])
				image[i] = sRGB1_color

		if dims == 2:
			for i in range(z_ar.shape[0]):
				for j in range(z_ar.shape[1]):
					a = self.a_func_of_theta(theta_ar[i,j])
					b = self.b_func_of_theta(theta_ar[i,j])

					lab_color = z_ar[i,j]*np.array([self.L_max,a,b])
					if scale_alpha: lab_color = np.array([self.L_max,a,b])

					sRGB1_color = cspace_convert(lab_color, self.color_space , "sRGB1")
					if clip_values:
						if any(sRGB1_color>1) or any(sRGB1_color<0):
							sRGB1_color_clipped = np.clip( sRGB1_color, 0,1)
							Lab_color_clipped = cspace_convert(sRGB1_color_clipped, "sRGB1", self.color_space )
							if verbose:print ('L*a*b*/sRGB1 color', np.array(lab_color), '/', sRGB1_color,
								'\n-----> clipped to',  Lab_color_clipped, '/', sRGB1_color_clipped)
							sRGB1_color = sRGB1_color_clipped

					if scale_alpha: sRGB1_color =  np.append(sRGB1_color, z_ar[i,j])
					image[i,j] = sRGB1_color

		return image


class isoluminant_uniform_circle_colormap:
	def __init__(self, radius = 40, L_max = 74, theta0 = 0.0, use_max_radius = True, color_space = default_colorspace):
		self.radius = radius
		self.L_max = L_max
		self.theta0 = theta0
		self.color_space = color_space
		### these are for compatilibilty
		self.a_knots = None
		self.b_knots = None
		self.min_theta_pts = None

		if use_max_radius:
			self.maximize_radius

	def maximize_radius(self,cut_off = 0.1, verbose=False):
		radius_max = 128.0
		radius_min = 0.0
		dE = 2*pi*cut_off # same spatial sampling, or something, I just made this up, it seems to work
		while (radius_max - radius_min) > cut_off:
			self.radius = 0.5*(radius_max + radius_min)

			circumference = self.radius * 2*pi
			ncirc = circumference/dE
			theta_levels =linspace(0,1,ncirc)

			sequence = self(theta_levels, np.ones(theta_levels.shape), clip_values = False)

			if np.any(sequence>1.0) or np.any(sequence<0.0):
				radius_is_safe = False
				if verbose:print (radius_min, self.radius, radius_max, radius_is_safe)
				radius_max = self.radius

			else:
				radius_is_safe = True
				if verbose:print (radius_min, self.radius, radius_max, radius_is_safe)
				radius_min = self.radius

		if radius_is_safe==False:
			self.radius = radius_min

	def a_func_of_theta(self, theta):
		return  self.radius*np.cos(theta+ self.theta0)
		#self.a_func_of_theta = a_func_of_theta

	def b_func_of_theta(self, theta):
		return  self.radius*np.sin( theta+self.theta0)

	def __call__(self, theta, z, scale_alpha = False, clip_values = True, verbose = True ):
		z_ar = np.array(z)
		theta_ar =  np.array(theta)

		dims = len(z_ar.shape)

		if scale_alpha:
			image = np.zeros(list(z_ar.shape)+[4])
		else:
			image = np.zeros(list(z_ar.shape)+[3])


		if dims == 1:
			for i in range(z_ar.shape[0]):
					a = self.a_func_of_theta(theta_ar[i] )
					b = self.b_func_of_theta(theta_ar[i] )

					lab_color = z_ar[i]*np.array([self.L_max,a,b])
					if scale_alpha: lab_color = np.array([self.L_max,a,b])

					sRGB1_color = cspace_convert(lab_color, self.color_space , "sRGB1")
					if clip_values:
						if any(sRGB1_color>1) or any(sRGB1_color<0):
							sRGB1_color_clipped = np.clip( sRGB1_color, 0,1)
							Lab_color_clipped = cspace_convert(sRGB1_color_clipped, "sRGB1", self.color_space )
							print (self.color_space +'/sRGB1 color', np.array(lab_color), '/', sRGB1_color,
								'\n-----> clipped to',  Lab_color_clipped, '/', sRGB1_color_clipped)
							sRGB1_color = sRGB1_color_clipped

					if scale_alpha: sRGB1_color =  np.append(sRGB1_color, z_ar[i])
					image[i] = sRGB1_color


		if dims == 2:
			for i in range(z_ar.shape[0]):
				for j in range(z_ar.shape[1]):
					a = self.a_func_of_theta(theta_ar[i,j] )
					b = self.b_func_of_theta(theta_ar[i,j] )

					lab_color = z_ar[i,j]*np.array([self.L_max,a,b])
					if scale_alpha: lab_color = np.array([self.L_max,a,b])

					sRGB1_color = cspace_convert(lab_color, self.color_space , "sRGB1")
					if clip_values:
						if any(sRGB1_color>1) or any(sRGB1_color<0):
							sRGB1_color_clipped = np.clip( sRGB1_color, 0,1)
							Lab_color_clipped = cspace_convert(sRGB1_color_clipped, "sRGB1", self.color_space )
							print ('L*a*b*/sRGB1 color', np.array(lab_color), '/', sRGB1_color,
								'\n-----> clipped to',  Lab_color_clipped, '/', sRGB1_color_clipped)
							sRGB1_color = sRGB1_color_clipped

					if scale_alpha: sRGB1_color =  np.append(sRGB1_color, z_ar[i,j])
					image[i,j] = sRGB1_color

		return image

def HSV_colormap(theta, z, scale_alpha = False, clip_values = True, verbose = True):
	from colorsys import hsv_to_rgb, rgb_to_hsv
	z_ar = np.array(z)
	theta_ar =  np.array(theta)/(2*np.pi)

	dims = len(z_ar.shape)

	if scale_alpha:
		image = np.zeros(list(z_ar.shape)+[4])
	else:
		image = np.zeros(list(z_ar.shape)+[3])


	if dims == 1:
		for i in range(z_ar.shape[0]):

			HSV_color = [theta_ar[i] ,1.0, z_ar[i] ]
			if scale_alpha: HSV_color[2] = 0.0
			sRGB1_color = np.array(hsv_to_rgb(*HSV_color))

			if clip_values:
				if any(sRGB1_color>1) or any(sRGB1_color<0):
					sRGB1_color_clipped = np.clip( sRGB1_color, 0,1)
					HSV_color_clipped = rgb_to_hsv(*sRGB1_color_clipped)
					if verbose:print ('HSV/sRGB1 color', HSV_color, '/', sRGB1_color,
							'\n-----> clipped to',  HSV_color_clipped, '/', sRGB1_color_clipped)
					sRGB1_color = sRGB1_color_clipped

			if scale_alpha: sRGB1_color =  np.append(sRGB1_color, z_ar[i])
			image[i] = sRGB1_color


	if dims == 2:
		for i in range(z_ar.shape[0]):
			for j in range(z_ar.shape[1]):

				HSV_color = [theta_ar[i,j] ,1.0, z_ar[i,j] ]
				if scale_alpha: HSV_color[2] = 0.0
				sRGB1_color = np.array(hsv_to_rgb(*HSV_color))

				if clip_values:
					if any(sRGB1_color>1) or any(sRGB1_color<0):
						sRGB1_color_clipped = np.clip( sRGB1_color, 0,1)
						HSV_color_clipped = rgb_to_hsv(*sRGB1_color_clipped)
						if verbose:print ('HSV/sRGB1 color', HSV_color, '/', sRGB1_color,
								'\n-----> clipped to',  HSV_color_clipped, '/', sRGB1_color_clipped)
						sRGB1_color = sRGB1_color_clipped

				if scale_alpha: sRGB1_color =  np.append(sRGB1_color, z_ar[i,j])
				image[i,j] = sRGB1_color


	return image

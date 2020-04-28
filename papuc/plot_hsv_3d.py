from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np

from papuc.core import default_colorspace

from colorspacious import cspace_convert

#####

# I'm not sure that 'antialiased' is needed to be off, but it's suggested online
# there needs to be *some* line width or the vectorized representations like pdf have stupid lines
# the joinstyle of bevel or round makes the lines not shoot outwards like miter
default_triangle_style = {'antialiased': False, 'linewidth' : 0.1, 'joinstyle': 'bevel' }


def hsv_color_triangle_mesh(
                num_theta_levels = 32,
                num_z_levels = 16,
                magnitude = 1.0,
                colorspace = default_colorspace,
                verbose = False):

    ''' color_space can be "CIELab" or "CAM02-UCS" '''

    from colorsys import hsv_to_rgb, rgb_to_hsv

    def eval_to_lab(theta, mag):
            HSV_color = [theta/(2*np.pi) ,1.0, mag ]
            sRGB1_color = np.array(hsv_to_rgb(*HSV_color))
            lab_color = cspace_convert(sRGB1_color, "sRGB1", colorspace)
            return lab_color

    def eval_vert_to_lab(vert):
        out = np.zeros((len(vert),3))
        for i in range(len(vert)):
            theta = vert[i][0]
            mag   = vert[i][1]
            HSV_color = [theta/(2*np.pi) ,1.0, mag ]
            sRGB1_color = np.array(hsv_to_rgb(*HSV_color))
            lab_color = cspace_convert(sRGB1_color, "sRGB1", colorspace)
            out[i] = lab_color
        return out


    mode = '2D'
    mode.upper()

    if mode.upper() == '2D':
        width = 2*np.pi
        height = magnitude
        delta = width/num_theta_levels
        deltay = -height/num_z_levels



    vert_list = []
    color_list = []

    if True:
        for iy in range(num_z_levels):
            for ix in range(num_theta_levels):
                # low
                theta = delta*(1/3.0+ix)
                mag   = height + deltay*(iy+1/3.0)
                color_list.append(eval_to_lab(theta,mag))
                vertices = [
                            (delta*(ix    ),  height + deltay*(iy)),
                            (delta*(1.0+ix),  height + deltay*(iy)),
                            (delta*(ix    ),  height + deltay*(iy+1.0))]
                vert_list.append(eval_vert_to_lab(vertices))

                # high
                theta = delta*(2/3.0+ix)
                mag   = height + deltay*(iy+2/3.0)
                color_list.append(eval_to_lab(theta,mag))
                vertices = [
                            (delta*(1.0+ix),  height + deltay*(iy    )),
                            (delta*(1.0+ix),  height + deltay*(iy+1.0)),
                            (delta*(ix    ),  height + deltay*(iy+1.0))]
                vert_list.append(eval_vert_to_lab(vertices))

    else: # this was a strange idea below
        if verbose: print('downwards triangles')
        for iy in range(nsegs):
            for ix in range(nsegs-iy):
                #f2 = (ix+0.5)*(1.0/nsegs)
                #f3 = (iy+0.5)*(1.0/nsegs)
                #fracs = np.array([1-f2-f3,f2,f3])
                #fracs = (1/nsegs) * np.array([nsegs-ix-1/3.-iy-1/3., ix+1/3., iy+1/3.])
                #color = fracs.dot(cp)

                theta = delta*(0.5+0.5*iy+ix)
                mag   = height + deltay*(iy+1/3.)



                if mode.upper() =='2D':
                    vertices = [
                            (delta*(    0.5*iy+ix),  height + deltay*   iy),
                            (delta*(0.5+0.5*iy+ix),  height + deltay*(iy+1)),
                            (delta*(1.0+0.5*iy+ix),  height + deltay*   iy)]
                else:
                    fracs = (1/nsegs) * np.array([
                                        [nsegs-ix-iy,   ix,   iy  ],
                                        [nsegs-ix-1-iy, ix+1, iy  ],
                                        [nsegs-ix-iy-1, ix  , iy+1]])
                    vertices = fracs.dot(cp)
                vert_list.append(eval_vert_to_lab(vertices))
                color_list.append(eval_to_lab(theta,mag))

        if verbose: print('upwards triangles')
        for iy in range(nsegs-1):
            for ix in range(nsegs-iy-1):
                #the centroids of 45-45-90 right triangles in this space
                #fracs = (1/nsegs) * np.array([nsegs-ix-2./3-iy-2./3, ix+2./3, iy+2./3])
                #color = fracs.dot(cp)

                theta = delta*(1.0+0.5*iy+ix)
                mag   = height + deltay*(iy+2/3.)

                if mode.upper()=='2D':
                    vertices = [
                                (delta*(0.5+0.5*iy+ix), height + deltay *(iy+1)),
                                (delta*(1.0+0.5*iy+ix), height + deltay * iy   ),
                                (delta*(1.5+0.5*iy+ix), height + deltay *(iy+1))]
                else:
                    fracs = (1/nsegs) * np.array([
                                        [nsegs-ix-1.0-iy,     ix+1.0,  iy    ],
                                        [nsegs-ix-iy-1.0,         ix,  iy+1.0],
                                        [nsegs-ix-1.0-iy-1.0, ix+1.0,  iy+1.0]])
                    vertices = fracs.dot(cp)

                #vert_list.append(vertices)
                #color_list.append(color)
                vert_list.append(eval_vert_to_lab(vertices))
                color_list.append(eval_to_lab(theta,mag))

    return vert_list, color_list




###########################
def hsv_UCS_cone_3D(ax,
                num_theta_levels = 32, num_z_levels = 16,
                #labels = ['Part 1','Part 2','Part 3' ],
                #order = [0,1,2],
                magnitude = 1.0,
                #font_options = {},
                #label_offset = 1.0,
                use_deuteranomaly = False, verbose = False,
                undisplayable_action = 'clip',
                color_for_undisplayable = (1.0, 1.0, 1.0, 0.0),
                triangle_style = default_triangle_style,
                colorspace = default_colorspace,
                alpha = 1.0 ):

    assert undisplayable_action in ['remove', 'clip', 'replace']
    #from matplotlib.patches import Polygon

    def fix_vert_order(verts):
        order = [1,2,0]
        verts_out  = np.zeros((3,3))
        for i in range(3):
            verts_out[i] = verts[i,order]
        return verts_out


    vert_list, color_list = hsv_color_triangle_mesh(
                num_theta_levels = num_theta_levels,
                num_z_levels = num_z_levels,
                magnitude = magnitude,
                colorspace = colorspace,
                verbose = verbose)

    #print(vert_list)
    #print(color_list)

    #print(norm)

    #print(color_list[-1])
    ## now we filter them
    vert_list_out = []
    color_list_out = []
    for vertices, color in zip(vert_list, color_list):
        sRGB1_color = cspace_convert(color, colorspace, "sRGB1")

        if use_deuteranomaly: sRGB1_color = cspace_convert(sRGB1_color, cvd_space, "sRGB1")
        if  np.any(sRGB1_color<0,0) or np.any(1.0<sRGB1_color):
            if undisplayable_action != 'remove':
                if undisplayable_action == 'replace':
                    sRGB1_color = color_for_undisplayable
                else:
                    sRGB1_color = np.clip(sRGB1_color,0,1)
                vert_list_out.append(fix_vert_order(vertices))
                color_list_out.append(sRGB1_color)

        else:
#            ax.add_patch(  Polygon(vertices, color =np.clip(sRGB1_color,0,1), **triangle_style))
            vert_list_out.append(fix_vert_order(vertices))
            color_list_out.append(sRGB1_color)

    poly = Poly3DCollection(vert_list_out, facecolors=color_list_out, edgecolors = color_list_out,  alpha=alpha, **triangle_style)
    ax.add_collection3d(poly, zdir='z')
    #art3d.Poly3DCollection(vert_list_out)

    #for vertices, sRGB1_color in zip(vert_list_out, color_list_out):
    #    ax.add_patch(  Polygon(vertices, color =sRGB1_color, **triangle_style))


    #if len(labels) ==3:
    #    cp = chemical_map.get_color_points()
    #    for i in range(3):
    #        ax.text(cp[i,1], cp[i,2], cp[i,0]+ label_offset,labels[i] , ha = 'center', va = 'bottom', **font_options)
        #ax.text(   0,             -label_offset, labels[0], ha = 'center', va = 'top',    **font_options)
        #ax.text(  width,         -label_offset, labels[1], ha = 'center', va = 'top',    **font_options)
        #ax.text( width/2.0, height+label_offset, labels[2], ha = 'center', va = 'bottom', **font_options)

    ax.set_xlabel('a*')
    ax.set_ylabel('b*')
    ax.set_zlabel('L*')



    L_test = np.array(vert_list)[:,:,0]
    a_test = np.array(vert_list)[:,:,1]
    #print(np.array(vert_list))
    #print(vert_list.shape)
    
    b_test = np.array(vert_list)[:,:,2]


    ax.set_xlim(a_test.min(), a_test.max())
    ax.set_ylim(b_test.min(), b_test.max())
    ax.set_zlim(L_test.min(), L_test.max())


    #ax.set_aspect('equal','box')
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,sin(60*pi/180))
    #ax.set_axis_off()

#############################

if __name__ == "__main__":

    from numpy import pi
    #from papuc.example_maps import colormaps
    #my_map = colormaps['default']


    #################
    fig = plt.figure()
    ax = fig.gca(projection='3d')


    hsv_UCS_cone_3D(ax, colorspace = "CIELab" )

    plt.show()

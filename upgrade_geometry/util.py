import numpy as np


def tilt_angle(nx,ny,nz):
    '''
    calculate tilt angle
    '''
    return np.arctan2(np.sqrt(nx*nx + ny*ny),nz)

def to_spherical(x,y,z):
    '''
    converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arctan2(np.sqrt(x*x + y*y),z)
    phi = np.arctan2(y,x)
    #####to ensure phi is positive
    if phi < 0.0:
        phi_pos = phi + 2*np.pi
    else:
        phi_pos = phi
    return r,theta,phi_pos


def to_spherical_lists(x_list,y_list,z_list):
    r_list = []
    θ_list = []
    φ_list = []
    for x,y,z in zip(x_list,y_list,z_list):
        r,theta,phi = to_spherical(x,y,z)
        r_list.append(r)
        θ_list.append(theta)
        φ_list.append(phi)
    return r_list,θ_list,φ_list
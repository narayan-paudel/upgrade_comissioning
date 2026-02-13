#!/usr/bin/env python

"""
Original file is located at
    https://colab.research.google.com/drive/1U5gYfEphIj7qy3DGPHenJ4kBLxQUjuN9
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

from pathlib import Path
home = str(Path.home())

plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"

#ellipse: AX2+BXY+CY2+DX+EY+F=0

def get_ellipse(h,k,a,b):
  x = []
  y = []
  # for iangle in np.linspace(0,2*np.pi,181):
  for iangle in np.linspace(0,2*np.pi,37):
    x.append(h + a * np.cos(iangle))
    y.append(k + b*np.sin(iangle))
  return np.asarray(x),np.asarray(y)

def rotate(x,y,theta):
  return x*np.cos(theta)-y*np.sin(theta),x*np.sin(theta)+y*np.cos(theta)


def get_rotated_ellipse(h,k,a,b,theta):
  x,y = get_ellipse(0,0,a,b)
  x,y = rotate(x,y,theta)
  return x+h,y+k

def get_noisy_ellipse(h,k,a,b,theta,noise_scale):
  x,y = get_rotated_ellipse(h,k,a,b,theta)
  x_noisy = [i+np.random.normal(0,noise_scale*abs(np.mean(x))) for i in x]
  y_noisy = [i+np.random.normal(0,noise_scale*abs(np.mean(y))) for i in y]
  return np.asarray(x_noisy), np.asarray(y_noisy)

x,y = get_noisy_ellipse(5,5,8,12,2*np.pi/7,0.03)


def fit_ellipse(x,y):
  x = np.asarray(x)
  y = np.asarray(y)
  A = np.stack([x**2,x*y,y**2,x,y]).T
  b = np.ones_like(x)
  w = np.linalg.lstsq(A, b, rcond=None)[0].squeeze()
  return w

def get_origin(A,B,C,D,E):
  return (2*C*D - B*E)/(B**2-4*A*C), (2*A*E-B*D)/(B**2-4*A*C)


def get_rotation(A,B,C,D,E):
  return 0.5*np.arctan2(-B,C-A)

def get_axes(A,B,C,D,E,F):
  t1 = B**2-4*A*C
  t2 = 2*(A*E**2+C*D**2-B*D*E+(B**2-4*A*C)*F)
  t3 = ((A+C)+np.sqrt((A-C)**2+B**2))
  t4 = ((A+C)-np.sqrt((A-C)**2+B**2))
  print(t1,t2,t3,t4)
  print(t2*t3/t1)
  return -np.sqrt(t2*t3)/t1, -np.sqrt(t2*t4)/t1

def get_eclipse_parameters(w):
  a,b = get_axes(*w,-1)
  theta = get_rotation(*w)
  x0,y0 = get_origin(*w)
  return np.asarray([x0,y0]),a,b,theta

def circle_from_ellipse(x,y,w):
  c0,a,b,theta = get_eclipse_parameters(w)
  r = np.sqrt(a*b)
#   ct, st = np.cos(theta), np.sin(theta)
  ct, st = np.cos(np.pi/2-theta), np.sin(np.pi/2-theta) #this works for multiple rotation
  R = np.array([[ct, -st],
                  [st,  ct]]) #rotation matrix for ellipse
  S = np.diag([a/r, b/r]) #scale for semi-major and semiminor axes
  A = (R.T)@S@R #to rotate the ellipse to align axes,scale to circle and rotate back
  x_circ = []
  y_circ = []
  for ix,iy in zip(x,y):
    ix_circ,iy_circ = (A@ (np.asarray([ix,iy])-c0).T).T
    x_circ.append(ix_circ)
    y_circ.append(iy_circ)
  return x_circ,y_circ

def corrected_ellipse(x,y):
  w = fit_ellipse(x,y)
  print(f"ellipse parameters:{w}")
  x_circ,y_circ = circle_from_ellipse(x,y,w)
  return x_circ,y_circ

def plot_corrected_circle(x,y):
  w = fit_ellipse(x,y)
  x_circ,y_circ = circle_from_ellipse(x,y,w)
  xlin = np.linspace(-10, 20, 30)
  ylin = np.linspace(-7, 17, 30)
  X, Y = np.meshgrid(xlin, ylin)

  Z = w[0]*X**2 + w[1]*X*Y + w[2]*Y**2 + w[3]*X + w[4]*Y

  fig, ax = plt.subplots()
  ax.scatter(x, y, label=f"ellipse data")
  ax.contour(X, Y, Z, [1],label=f"fit")
  ax.scatter(x_circ, y_circ, label=f"corrected circle")
  ax.set_aspect("equal")
  ax.grid()
  ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
  ax.grid(True,alpha=0.6)
  ax.legend(ncols=1,fontsize=8)
  ax.set_xlabel(r" B$_{x}$ [$\mu$T]", fontsize=16)
  ax.set_ylabel(r" B$_{y}$ [$\mu$T]", fontsize=16)
  plt.legend()
  plt.savefig(plotFolder+f"/corrected_ellipse.png",transparent=False,bbox_inches='tight')
  plt.savefig(plotFolder+f"/corrected_ellipse.pdf",transparent=False,bbox_inches='tight')
  plt.close()

plot_corrected_circle(x,y)


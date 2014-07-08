# Filename: halfarDome.py
# This is a module used by both halfar.py and halfar_results.py for calculating 
# the analytic solution for the Halfar dome.
# It is not meant to be run on its own.

# Written by Matt Hoffman, October 2013

import numpy as np


# Define the function to calculate the Halfar thickness
# Halfar, P. 1983. On the Dynamics of the Ice Sheets 2.  Journal of Geophysical Research, 88, 6043-6051.
def halfarDome(t,x,y,flwa,rhoi):
  # Input: t - time in years
  # Input: x - 1-d array of cell center x-positions
  # Input: y - 1-d array of cell center y-positions
  # Input: flwa - flow law parameter A in units of Pa^-3 yr^-1
  # Input: rhoi - ice density in kg/m3

  # Initial radius and central thickness of dome
  R0 = 60000.0 * np.sqrt(0.125)
  H0 = 2000.0 * np.sqrt(0.125)

  n = 3.0
  grav = 9.8101
  alpha = 1.0/9.0
  beta = 1.0/18.0
  secpera = 31556926.0
  Gamma = 2.0/(n+2.0) * flwa * (rhoi * grav)**n

  xcenter = max(x)/2.0
  ycenter = max(y)/2.0

  t0 = (beta/Gamma) * (7.0/4.0)**3 * (R0**4/H0**7)  # Note: this line assumes n=3!
  tr=(t+t0)/t0 

  H=np.zeros((len(y), len(x)))
  for i in range(len(x)):
    for j in range(len(y)):
      r = np.sqrt( (x[i]-xcenter)**2 + (y[j]-ycenter)**2)
      r=r/R0
      inside = max(0.0, 1.0 - (r / tr**beta)**((n+1.0) / n))

      H[j,i] = H0 * inside**(n / (2.0*n+1.0)) / tr**alpha
  return H.astype(np.float32)

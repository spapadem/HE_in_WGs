import numpy as np
from freq_used import *
## Problem setup.
# Wavenumber & source.
c0 = 1500. # Constant wave speed.
f0 = 75.   # Reference frequency.
lambda_0 = c0/f0; # Reference wavelength. We use this to define all sizes with respect to this scale.


# Waveguide characteristics
Dc = 7*lambda_0 # Waveguide constant depth
Dm = 20*lambda_0 # Waveguide max depth
Wc =  10*lambda_0 # Width with constant depth
Wm = 35*lambda_0 # Waveguide max width. (originally infinite but we have to truncate for computational purposes).
PML_size = 4*lambda_0 # Length of the Perfectly Matched Layer (PML) that helps us truncate our computational domain.

# Location and size of the scatterer.
x_sc = 19.5*lambda_0 # Location in x-axis
y_sc = 5*lambda_0 # Location in y-axis
b = 1*lambda_0 # Size of the scatterer (radius)

# # Source present in the waveguide.
# frq = 33. # Frequency in which the source emits its pulse.
omega = 2.*np.pi*frq # Angular frequency.
k = omega/c0 # wavenumber
x_s= Wm-3*lambda_0 # Position of source in x-axis.
Nr = 41
h = Dm/(Nr-1)
y_a = np.linspace(h, Dm-h, Nr)
r = 5 # Radius of source
alpha = np.log(10^6)/r**2


Nx = 201 # Number of points in the regular grid, in the x-direction.
Ny = 201 # Number of points in the regular grid, in the y-direction.
xg = np.linspace(x_sc-4*lambda_0,x_sc+4*lambda_0,Nx)
yg = np.linspace(y_sc-4*lambda_0,y_sc+4*lambda_0,Ny)

fmin = 41
fmax = 89
h = 2
Nf = int((fmax-fmin)/h + 1)
frequencies = np.linspace(fmin, fmax, Nf)
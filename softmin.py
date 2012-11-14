##################################################
## Script for minimising soft sphere potential  ## 
## for polydisperse disks in 2D.                ##
## Copyright 2012 D. Asenjo & F. Paillusson.    ##
##################################################

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#                            
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from softmin.soft_sphere import SoftSphere

npart = 64     # Number of particles.
phi = 0.9      # Volume fraction.
poly = 0.2     # Polydispersity.

mean = np.sqrt(phi/(np.pi*npart*(poly*poly+1.0)))    # Mean radius.
diams = np.random.normal(mean,poly*mean,[npart])     # Pick radii from Gaussian distribution. 
diams *= 2.0 

print 'Real volume fraction:', np.sum(0.25*np.pi*diams**2)

pot = SoftSphere(diams = diams)                      # Use soft sphere harmonic potential.

coords = np.random.uniform(-0.5,0.5,[npart*2])       # Start from random positions.
 
nsteps = 100000    # Maximum number of steps for steepest descent.
tol = 1.0e-7       # Stop when gradient < tol.
dx = 1.0e-3        # Stepsize for steepest descent

E,grad = pot.getEnergyGradient(coords)     # Calculate potential energy and gradient.

# Steepest descent loop.
for k in range(nsteps):            
    
    coords -= grad*dx
    E,grad = pot.getEnergyGradient(coords)
    rms = np.linalg.norm(grad)/np.sqrt(len(grad))
    if rms < tol:
        break

# Print final coords and diameters to file.
newcoords = np.reshape(coords,(npart,2))
np.savetxt('coords.dat',newcoords)
np.savetxt('diams.dat',diams)
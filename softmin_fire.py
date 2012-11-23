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

dt = 0.1
maxmove = 0.1
dtmax = 1.0
Nmin = 5
finc = 1.1
fdec = 0.5
astart = 0.1
fa = 0.99
a = 0.1
steps = 0
Nsteps = 0
v = np.zeros(2*npart)

E,grad = pot.getEnergyGradient(coords)     # Calculate potential energy and gradient.

# FIRE loop.
for k in range(nsteps):            
    
    E,grad = pot.getEnergyGradient(coords)
    f = np.sqrt(np.vdot(grad,grad))
    
    if f < tol:
        break
    
    P = np.vdot(-1.0*grad,v)

    if P > 0.0:
        v = (1.0-a)*v + a*(-1.0*grad/f)*np.sqrt(np.vdot(v,v))
        
        if (Nsteps > Nmin):
            dt = min(dt*finc,dtmax)
            a *= fa
        Nsteps = Nsteps + 1
    else:

        v = np.zeros(2*npart)
        a = astart
        dt = dt*fdec
        Nsteps = 0

    v = v - dt*grad
    dr = dt*v
    normdr = np.sqrt(np.vdot(dr,dr))

    if normdr > maxmove: 
        dr = dr * maxmove/normdr

    coords += dr

# Print final coords and diameters to file.
newcoords = np.reshape(coords,(npart,2))
np.savetxt('coords.dat',newcoords)
np.savetxt('diams.dat',diams)

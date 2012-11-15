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
from softmin.soft_sphere_pot import soft_sphere_pot

__all__ = ["SoftSphere"]

class SoftSphere():
    def __init__(self, diams, dimen = 2 ):
        self.diams = diams
        self.dimen = dimen
        
        dmean = np.mean(diams)
        dstd = np.std(diams)
        print "using soft sphere potential with mean diameter", dmean, "st. dev.", dstd
        
    
    def getEnergy(self, coords):
        natoms = len(coords)/self.dimen
        energy, force = soft_sphere_pot(self.dimen, coords, self.diams, [natoms])
        return energy
    
    def getEnergyGradient(self, coords):
        natoms = len(coords)/self.dimen
        energy, force = soft_sphere_pot(self.dimen, coords, self.diams, [natoms])
        return energy, force


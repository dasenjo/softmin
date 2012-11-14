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


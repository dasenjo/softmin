from numpy.distutils.core import setup
from numpy.distutils.core import Extension

setup(
packages=["softmin"], 
ext_modules=[ Extension("softmin.soft_sphere_pot", ["softmin/soft_sphere_pot.f90"])]
)

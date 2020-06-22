"""
To compile: 

1. make sure that setup.py reads in the current version (here simlib.pyx) 
2. python setup.py build_ext --inplace

"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[ Extension("simlib",
              ["simlib.pyx"],
              libraries=["m"],
              extra_compile_args = ["-ffast-math"])]



setup(
  name = "simlib",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules)

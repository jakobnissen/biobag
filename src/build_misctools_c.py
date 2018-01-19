# python build_misctools_c.py build_ext --inplace
from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'misctools_c',
  ext_modules = cythonize("misctools_c.pyx"),
)

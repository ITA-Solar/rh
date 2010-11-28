from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
 
NUMPY_INCLUDE = '/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/'

# The ugly: libraries going as extra link args (I can't be bothered now to find out
# how to do this properly).
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("rhpy", ["rhpy.pyx","tiago_libs.c"],
                             extra_link_args=["../librh.a","../librh_f90.a","../rhsc2d/project.o"],
                             include_dirs=[NUMPY_INCLUDE])])

#!/bin/bash
echo 'gcc -c tiago.c -o tiago.o'
gcc -c tiago_libs.c -o tiago_libs.o
echo 'cython rhpy.pyx'
cython rhpy.pyx
echo 'gcc -c -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/ rhpy.c -o rhpy.o'
gcc -c -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/ rhpy.c -o rhpy.o
echo 'gcc -bundle -undefined dynamic_lookup rhpy.o tiago.o rhsc2d/project.o librh.a librh_f90.a -o rhpy.so'
gcc -bundle -undefined dynamic_lookup rhpy.o tiago_libs.o ../rhsc2d/project.o ../librh.a ../librh_f90.a -o rhpy.so

#! /bin/bash
g++ -o SALTSA_example -O3 -funroll-loops -march=native -ftree-vectorize -Wall -fmessage-length=0 -fopenmp -fPIC -MMD -MP SALTSA_example.cpp SALTSA_astro.cpp SALTSA_file_functions.cpp SALTSA_interpolator.cpp SALTSA_orbit.cpp tk_spline.cpp

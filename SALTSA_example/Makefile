CXX = g++
CXXFLAGS = -O3 -Wall -O3 -funroll-loops -march=native -ftree-vectorize -Wall -fmessage-length=0 -fopenmp -fPIC -MMD -MP
LDFLAGS = -fopenmp

%.o:           	%.cpp
		$(CXX) $(CXXFLAGS) -c $*.cpp

everything:    	SALTSA_example 

SALTSA_example_obj = SALTSA_example.o SALTSA_interpolator.o tk_spline.o SALTSA_astro.o SALTSA_orbit.o SALTSA_file_functions.o

SALTSA_example:	$(SALTSA_example_obj)
		$(CXX) -o $@ $(SALTSA_example_obj) $(LDFLAGS) -L. -lm

SALTSA_example.o:	SALTSA_example.cpp SALTSA.h SALTSA_calculus.hpp SALTSA_file_functions.h SALTSA_global.h SALTSA_interpolator.h SALTSA_astro.h SALTSA_orbit.h tk_spline.h SALTSA_unitconvs.h SALTSA_cache.hpp SALTSA_functor.hpp SALTSA_misc_functions.hpp SALTSA_phase.hpp

SALTSA_interpolator.o:	SALTSA_interpolator.cpp SALTSA_interpolator.h tk_spline.h

tk_spline.o:   	tk_spline.cpp tk_spline.h

SALTSA_astro.o:	SALTSA_astro.cpp SALTSA_global.h SALTSA_astro.h SALTSA_misc_functions.hpp SALTSA_calculus.hpp SALTSA_solvers.hpp SALTSA_unitconvs.h SALTSA_cache.hpp SALTSA_functor.hpp SALTSA_file_functions.h SALTSA_vector.hpp

SALTSA_orbit.o:	SALTSA_orbit.cpp SALTSA_global.h SALTSA_unitconvs.h SALTSA_misc_functions.hpp SALTSA_interpolator.h SALTSA_astro.h SALTSA_file_functions.h SALTSA_phase.hpp SALTSA_functor.hpp SALTSA_calculus.hpp SALTSA_solvers.hpp SALTSA_orbit.h SALTSA.h tk_spline.h SALTSA_cache.hpp SALTSA_vector.hpp

SALTSA_file_functions.o:	SALTSA_file_functions.cpp SALTSA_misc_functions.hpp SALTSA_file_functions.h SALTSA_global.h

clean:
	rm -f *.o *.d

distclean:
	rm -f SALTSA_example *.o *.d

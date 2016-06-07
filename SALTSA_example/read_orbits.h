#ifndef INC_READ_ORBITS
#define INC_READ_ORBITS

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <list>

#include "orbit_class.h"
#include "interloper_class.h"

class Interloper;

#ifdef ORB_LIST
void read_orbit_data(const char*, std::list <Orbit>&);
#else
void read_orbit_data(const char*, std::vector<Orbit>&, unsigned int N_start=0, unsigned int N_to_read=0);
#endif

void read_interloper_data(const char*, std::vector<Interloper>&);

#endif //INC_READ_ORBITS

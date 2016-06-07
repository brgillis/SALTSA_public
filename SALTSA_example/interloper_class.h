#ifndef INC_INTERLOPER_CLASS
#define INC_INTERLOPER_CLASS

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "read_orbits.h"

class Interloper {
private:
	long int ID;
	float mass;
	float host_mass;
	float host_rvir;
	float host_vrms;
	float R;
	float V;
public:
	Interloper();
	Interloper(std::ifstream&);
	Interloper(const Interloper&);
	void swap(Interloper&);
#ifndef SWIG
	Interloper& operator=(Interloper);
#endif
	~Interloper();
	long int get_ID();
	float get_mass();
	float get_host_mass();
	float get_host_rvir();
	float get_host_vrms();
	float get_R();
	float get_V();
	void set_ID(long int);
	void set_mass(float);
	void set_host_mass(float);
	void set_host_rvir(float);
	void set_host_vrms(float);
	void set_R(float);
	void set_V(float);
};

#endif //INC_INTERLOPER CLASS

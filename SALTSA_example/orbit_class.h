#ifndef INC_ORBIT_CLASS
#define INC_ORBIT_CLASS

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <list>
#include "brg/physics/SALTSA/stripping_orbit.h"
#include "brg/physics/units/unit_conversions.hpp"

class Orbit {
private:
	//member variables
	int nsteps; //number of timesteps for this orbit
	std::vector<long int> id;
	std::vector<float> scale;
	std::vector<float> mass;
	std::vector<int> parent;
	std::vector<int> subsub;
	std::vector<float> parent_mass;
	std::vector<float> parent_radius;
	std::vector<float> snp_mass;
	std::vector<float> snp_radius;
	std::vector<float> snp_vrms;
	std::vector<int> fppb;
	std::vector<float> fppb_mass;
	std::vector<float> fppb_radius;
	std::vector<float> fppb_vrms;
	std::vector<int> fppb_eq_snp;
	std::vector<int> host_id;
	std::vector<float> xp;
	std::vector<float> yp;
	std::vector<float> zp;
	std::vector<float> vxp;
	std::vector<float> vyp;
	std::vector<float> vzp;
	std::vector<float> xf;
	std::vector<float> yf;
	std::vector<float> zf;
	std::vector<float> vxf;
	std::vector<float> vyf;
	std::vector<float> vzf;

	#if (1) //brg
	brgastro::stripping_orbit *orbit_spline_ptr;
	bool orbit_spline_loaded;
	#endif // brg

public:
	//constructors
	Orbit();
	Orbit(std::ifstream&);
	//copy constructor
	Orbit(const Orbit&);
	//assignment operator & swap
	void swap(Orbit&);
	Orbit& operator=(Orbit);
	//destructor
	~Orbit();
	//gets and sets
	int get_nsteps();
	long int get_id(int);
	float get_scale(int);
	float get_mass(int);
	bool get_parent(int);
	bool get_subsub(int);
	float get_parent_mass(int);
	float get_parent_radius(int);
	float get_snp_mass(int);
	float get_snp_radius(int);
	float get_snp_vrms(int);
	bool get_fppb(int);
	float get_fppb_mass(int);
	float get_fppb_radius(int);
	float get_fppb_vrms(int);
	bool get_fppb_eq_snp(int);
	int get_host_id(int);
	float get_xp(int);
	float get_yp(int);
	float get_zp(int);
	float get_vxp(int);
	float get_vyp(int);
	float get_vzp(int);
	float get_xf(int);
	float get_yf(int);
	float get_zf(int);
	float get_vxf(int);
	float get_vyf(int);
	float get_vzf(int);

	void set_nsteps(int);
	void set_id(int, long int);
	void set_scale(int, float);
	void set_mass(int, float);
	void set_parent(int, bool);
	void set_subsub(int, bool);
	void set_parent_mass(int, float);
	void set_parent_radius(int, float);
	void set_snp_mass(int, float);
	void set_snp_radius(int, float);
	void set_snp_vrms(int, float);
	void set_fppb(int, bool);
	void set_fppb_mass(int, float);
	void set_fppb_radius(int, float);
	void set_fppb_vrms(int, float);
	void set_fppb_eq_snp(int, bool);
	void set_host_id(int, int);
	void set_xp(int, float);
	void set_yp(int, float);
	void set_zp(int, float);
	void set_vxp(int, float);
	void set_vyp(int, float);
	void set_vzp(int, float);
	void set_xf(int, float);
	void set_yf(int, float);
	void set_zf(int, float);
	void set_vxf(int, float);
	void set_vyf(int, float);
	void set_vzf(int, float);

	//other functions
	void display();
	float rp(int);
	float vp(int);
	float rf(int);
	float vf(int);
	float vrp(int);
	float vrf(int);
	float vtp(int);
	float vtf(int);
	float rpproj(int, int);
	float vpproj(int, int);
	float rfproj(int, int);
	float vfproj(int, int);

	#if (1) //brg
	brgastro::stripping_orbit *stripping_orbit_spline();
	void unload_orbit_spline();
	const int load_orbit_spline(const int init_step, brgastro::stripping_orbit *new_ptr, const bool silent=false);
	#endif // brg
};

#ifdef SWIG
#ifdef ORB_LIST
Orbit extract_orbit(std::list <Orbit>&, int);
#else
Orbit extract_orbit(std::vector <Orbit>&, int);
#endif
#ifdef ORB_LIST
int get_size(std::list <Orbit>&);
#else
int get_size(std::vector <Orbit>&);
#endif
#endif

#endif //INC_ORBIT_CLASS

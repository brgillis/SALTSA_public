#include "orbit_class.h"
#include "brg/physics/astro.h"

using namespace std;

namespace unitconv = brgastro::unitconv;

#ifdef SWIG
//global var for use in python interface
#ifdef ORB_LIST
list<Orbit> py_orbits = list<Orbit>();
#else
vector<Orbit> py_orbits = vector<Orbit>();
#endif

//function to pull Orbit objects from a list<Orbit> for use in python interface
#ifdef ORB_LIST
Orbit extract_orbit(list <Orbit>& Orbits, int n)
{
	int count = 0;
	list<Orbit>::iterator itty;
	for(itty = Orbits.begin(); itty != Orbits.end(); itty++, count++)
	{
		if(count == n) return *itty;
	}
	cout << "Failed to return Orbit from list: end of list reached." << endl;
	exit(0);
}
#else
Orbit extract_orbit(vector <Orbit>& Orbits, int n)
{
	return Orbits[n];
}
#endif

#ifdef ORB_LIST
int get_size(list <Orbit>& Orbits)
{
	return Orbits.size();
}
#else
int get_size(vector <Orbit>& Orbits)
{
	return Orbits.size();
}
#endif

#endif

//constructors
Orbit::Orbit() {
	nsteps = 0;
	id.clear();
	scale.clear();
	mass.clear();
	parent.clear();
	subsub.clear();
	parent_mass.clear();
	parent_radius.clear();
	snp_mass.clear();
	snp_radius.clear();
	snp_vrms.clear();
	fppb.clear();
	fppb_mass.clear();
	fppb_radius.clear();
	fppb_vrms.clear();
	fppb_eq_snp.clear();
	host_id.clear();
	xp.clear();
	yp.clear();
	zp.clear();
	vxp.clear();
	vyp.clear();
	vzp.clear();
	xf.clear();
	yf.clear();
	zf.clear();
	vxf.clear();
	vyf.clear();
	vzf.clear();

	#if (1) //brg
	orbit_spline_loaded = false;
	orbit_spline_ptr = 0;
	#endif // brg
	return;
}

Orbit::Orbit(ifstream& data) {
	#if (1) //brg
	orbit_spline_loaded = false;
	orbit_spline_ptr = 0;
	#endif // brg
	string junk;
	float buffer;
	long int ibuffer;
	int ns;
	int nf;
	data >> ns;
	data >> nf;
	if (ns >= nf)
		nsteps = ns;
	else
		nsteps = nf;
	id.reserve(nsteps);
	scale.reserve(nsteps);
	mass.reserve(nsteps);
	parent.reserve(nsteps);
	subsub.reserve(nsteps);
	parent_mass.reserve(nsteps);
	parent_radius.reserve(nsteps);
	snp_mass.reserve(nsteps);
	snp_radius.reserve(nsteps);
	snp_vrms.reserve(nsteps);
	fppb.reserve(nsteps);
	fppb_mass.reserve(nsteps);
	fppb_radius.reserve(nsteps);
	fppb_vrms.reserve(nsteps);
	fppb_eq_snp.reserve(nsteps);
	host_id.reserve(nsteps);
	xp.reserve(nsteps);
	yp.reserve(nsteps);
	zp.reserve(nsteps);
	vxp.reserve(nsteps);
	vyp.reserve(nsteps);
	vzp.reserve(nsteps);
	xf.reserve(nsteps);
	yf.reserve(nsteps);
	zf.reserve(nsteps);
	vxf.reserve(nsteps);
	vyf.reserve(nsteps);
	vzf.reserve(nsteps);
	for (int loop = 0; loop < nsteps; loop++) {
		//id
		data >> ibuffer;
		id.push_back(ibuffer);
		//scale
		data >> buffer;
		scale.push_back(buffer);
		//halo mass
		data >> buffer;
		mass.push_back(buffer);
		//parent mass
		data >> buffer;
		if (buffer <= 0) {
			parent.push_back(false);
			parent_mass.push_back(-1);
		} else {
			parent.push_back(true);
			parent_mass.push_back(buffer);
		}
		//parent radius
		data >> buffer;
		if (buffer <= 0)
			parent_radius.push_back(-1);
		else
			parent_radius.push_back(buffer);
		//subsub
		data >> buffer;
		if (buffer <= 0)
			subsub.push_back(false);
		else
			subsub.push_back(true);
		//snp mass
		data >> buffer;
		if (buffer <= 0)
			snp_mass.push_back(-1);
		else
			snp_mass.push_back(buffer);
		//snp radius
		data >> buffer;
		if (buffer <= 0)
			snp_radius.push_back(-1);
		else
			snp_radius.push_back(buffer);
		//snp vrms
		data >> buffer;
		if (buffer <= 0)
			snp_vrms.push_back(-1);
		else
			snp_vrms.push_back(buffer);
		//fppb mass
		data >> buffer;
		if (buffer <= 0) {
			fppb_mass.push_back(-1);
			fppb.push_back(false);
		} else {
			fppb_mass.push_back(buffer);
			fppb.push_back(true);
		}
		//fppb radius
		data >> buffer;
		if (buffer <= 0)
			fppb_radius.push_back(-1);
		else
			fppb_radius.push_back(buffer);
		//fppb vrms
		data >> buffer;
		if (buffer <= 0)
			fppb_vrms.push_back(-1);
		else
			fppb_vrms.push_back(buffer);
		//fppb_eq_snp
		data >> buffer;
		if (buffer <= 0)
			fppb_eq_snp.push_back(false);
		else
			fppb_eq_snp.push_back(true);
		//fppb_eq_snp
		data >> buffer;
		host_id.push_back(buffer);
		getline(data, junk);
		if (data.peek() == 'V') {
			getline(data, junk);
			xp.push_back(-1);
			yp.push_back(-1);
			zp.push_back(-1);
			vxp.push_back(-1);
			vyp.push_back(-1);
			vzp.push_back(-1);
		} else {
			//coordinates relative to parent
			data >> buffer;
			xp.push_back(buffer);
			data >> buffer;
			yp.push_back(buffer);
			data >> buffer;
			zp.push_back(buffer);
			data >> buffer;
			vxp.push_back(buffer);
			data >> buffer;
			vyp.push_back(buffer);
			data >> buffer;
			vzp.push_back(buffer);
			getline(data, junk);
		}
		if (data.peek() == 'V') {
			getline(data, junk);
			xf.push_back(-1);
			yf.push_back(-1);
			zf.push_back(-1);
			vxf.push_back(-1);
			vyf.push_back(-1);
			vzf.push_back(-1);
		} else {
			//coordinates relative to FPPB
			data >> buffer;
			xf.push_back(buffer);
			data >> buffer;
			yf.push_back(buffer);
			data >> buffer;
			zf.push_back(buffer);
			data >> buffer;
			vxf.push_back(buffer);
			data >> buffer;
			vyf.push_back(buffer);
			data >> buffer;
			vzf.push_back(buffer);
			getline(data, junk);
		}
	}
	return;
}

//copy constructor
Orbit::Orbit(const Orbit& other) {
	nsteps = other.nsteps;
	id = other.id;
	scale = other.scale;
	mass = other.mass;
	parent = other.parent;
	parent_mass = other.parent_mass;
	parent_radius = other.parent_radius;
	subsub = other.subsub;
	snp_mass = other.snp_mass;
	snp_radius = other.snp_radius;
	snp_vrms = other.snp_vrms;
	fppb = other.fppb;
	fppb_mass = other.fppb_mass;
	fppb_radius = other.fppb_radius;
	fppb_vrms = other.fppb_vrms;
	fppb_eq_snp = other.fppb_eq_snp;
	host_id = other.host_id;
	xp = other.xp;
	yp = other.yp;
	zp = other.zp;
	vxp = other.vxp;
	vyp = other.vyp;
	vzp = other.vzp;
	xf = other.xf;
	yf = other.yf;
	zf = other.zf;
	vxf = other.vxf;
	vyf = other.vyf;
	vzf = other.vzf;
	orbit_spline_loaded = other.orbit_spline_loaded;
	if (orbit_spline_loaded) {
		orbit_spline_ptr = other.orbit_spline_ptr->stripping_orbit_clone();
	} else {
		orbit_spline_ptr = 0;
	}

	return;
}

//swap & assignment operator
void Orbit::swap(Orbit& other) {
	std::swap(nsteps, other.nsteps);
	std::swap(id, other.id);
	std::swap(scale, other.scale);
	std::swap(mass, other.mass);
	std::swap(parent, other.parent);
	std::swap(parent_mass, other.parent_mass);
	std::swap(parent_radius, other.parent_radius);
	std::swap(subsub, other.subsub);
	std::swap(snp_mass, other.snp_mass);
	std::swap(snp_radius, other.snp_radius);
	std::swap(snp_vrms, other.snp_vrms);
	std::swap(fppb, other.fppb);
	std::swap(fppb_mass, other.fppb_mass);
	std::swap(fppb_radius, other.fppb_radius);
	std::swap(fppb_vrms, other.fppb_vrms);
	std::swap(fppb_eq_snp, other.fppb_eq_snp);
	std::swap(host_id, other.host_id);
	std::swap(xp, other.xp);
	std::swap(yp, other.yp);
	std::swap(zp, other.zp);
	std::swap(vxp, other.vxp);
	std::swap(vyp, other.vyp);
	std::swap(vzp, other.vzp);
	std::swap(xf, other.xf);
	std::swap(yf, other.yf);
	std::swap(zf, other.zf);
	std::swap(vxf, other.vxf);
	std::swap(vyf, other.vyf);
	std::swap(vzf, other.vzf);
	std::swap(orbit_spline_ptr, other.orbit_spline_ptr);
	std::swap(orbit_spline_loaded, other.orbit_spline_loaded);

	return;
}

Orbit& Orbit::operator=(Orbit other) {
	Orbit::swap(other);

	return *this;
}

//destructor
Orbit::~Orbit() {
	nsteps = 0;
	id.clear();
	scale.clear();
	mass.clear();
	parent.clear();
	parent_mass.clear();
	parent_radius.clear();
	subsub.clear();
	snp_mass.clear();
	snp_radius.clear();
	snp_vrms.clear();
	fppb.clear();
	fppb_mass.clear();
	fppb_radius.clear();
	fppb_vrms.clear();
	fppb_eq_snp.clear();
	host_id.clear();
	xp.clear();
	yp.clear();
	zp.clear();
	vxp.clear();
	vyp.clear();
	vzp.clear();
	xf.clear();
	yf.clear();
	zf.clear();
	vxf.clear();
	vyf.clear();
	vzf.clear();

	unload_orbit_spline();

	return;
}

//gets and sets
int Orbit::get_nsteps() {
	return nsteps;
}

long int Orbit::get_id(int step) {
	return id[step];
}

float Orbit::get_scale(int step) {
	return scale[step];
}

float Orbit::get_mass(int step) {
	return mass[step];
}

bool Orbit::get_parent(int step) {
	return bool(parent[step]);
}

bool Orbit::get_subsub(int step) {
	if (subsub[step] <= 0)
		return false;
	else
		return true;
}

float Orbit::get_parent_mass(int step) {
	return parent_mass[step];
}

float Orbit::get_parent_radius(int step) {
	return parent_radius[step];
}

float Orbit::get_snp_mass(int step) {
	return snp_mass[step];
}

float Orbit::get_snp_radius(int step) {
	return snp_radius[step];
}

float Orbit::get_snp_vrms(int step) {
	return snp_vrms[step];
}

bool Orbit::get_fppb(int step) {
	return fppb[step];
}

float Orbit::get_fppb_mass(int step) {
	return fppb_mass[step];
}

float Orbit::get_fppb_radius(int step) {
	return fppb_radius[step];
}

float Orbit::get_fppb_vrms(int step) {
	return fppb_vrms[step];
}

bool Orbit::get_fppb_eq_snp(int step) {
	return fppb_eq_snp[step];
}

int Orbit::get_host_id(int step) {
	return host_id[step];
}

float Orbit::get_xp(int step) {
	return xp[step];
}

float Orbit::get_yp(int step) {
	return yp[step];
}

float Orbit::get_zp(int step) {
	return zp[step];
}

float Orbit::get_vxp(int step) {
	return vxp[step];
}

float Orbit::get_vyp(int step) {
	return vyp[step];
}

float Orbit::get_vzp(int step) {
	return vzp[step];
}

float Orbit::get_xf(int step) {
	return xf[step];
}

float Orbit::get_yf(int step) {
	return yf[step];
}

float Orbit::get_zf(int step) {
	return zf[step];
}

float Orbit::get_vxf(int step) {
	return vxf[step];
}

float Orbit::get_vyf(int step) {
	return vyf[step];
}

float Orbit::get_vzf(int step) {
	return vzf[step];
}

void Orbit::set_nsteps(int n) {
	nsteps = n;
	return;
}

void Orbit::set_id(int step, long int i) {
	id[step] = i;
	return;
}

void Orbit::set_scale(int step, float a) {
	scale[step] = a;
	return;
}

void Orbit::set_mass(int step, float m) {
	mass[step] = m;
	return;
}

void Orbit::set_parent(int step, bool p) {
	parent[step] = p;
	if (p == true && (parent_mass[step] <= 0 || parent_radius[step] <= 0))
		cout
				<< "WARNING: Parent flag true with unphysical parent mass or radius."
				<< endl;
	if (p == false && (parent_mass[step] > 0 || parent_radius[step] > 0))
		cout
				<< "WARNING: Parent flag false with non-zero parent mass or radius."
				<< endl;
	return;
}

void Orbit::set_subsub(int step, bool ss) {
	subsub[step] = ss;
	if (ss == true && parent[step] == false)
		cout << "WARNING: Subsub flag true with parent flag false." << endl;
	if (ss == false && parent[step] == false)
		subsub[step] = -1;
	return;
}

void Orbit::set_parent_mass(int step, float pm) {
	parent_mass[step] = pm;
	if (pm > 0)
		parent[step] = 1;
	if (pm <= 0)
		parent[step] = 0;
	return;
}

void Orbit::set_parent_radius(int step, float pr) {
	parent_radius[step] = pr;
	if (pr > 0)
		parent[step] = 1;
	if (pr <= 0)
		parent[step] = 0;
	return;
}

void Orbit::set_snp_mass(int step, float snpm) {
	snp_mass[step] = snpm;
	return;
}

void Orbit::set_snp_radius(int step, float snpr) {
	snp_radius[step] = snpr;
	return;
}

void Orbit::set_snp_vrms(int step, float snpv) {
	snp_vrms[step] = snpv;
	return;
}

void Orbit::set_fppb(int step, bool f) {
	fppb[step] = f;
	if (f == true
			&& (fppb_mass[step] <= 0 || fppb_radius[step] <= 0
					|| fppb_vrms[step] <= 0))
		cout
				<< "WARNING: FPPB flag true with unphysical FPPB mass or radius or vrms."
				<< endl;
	if (f == false
			&& (fppb_mass[step] > 0 || fppb_radius[step] > 0
					|| fppb_vrms[step] > 0))
		cout
				<< "WARNING: FPPB flag false with non-zero FPPB mass or radius or vrms."
				<< endl;
	return;
}

void Orbit::set_fppb_mass(int step, float fm) {
	fppb_mass[step] = fm;
	if (fm > 0)
		fppb[step] = 1;
	if (fm <= 0)
		fppb[step] = 0;
	return;
}

void Orbit::set_fppb_radius(int step, float fr) {
	fppb_radius[step] = fr;
	if (fr > 0)
		fppb[step] = 1;
	if (fr <= 0)
		fppb[step] = 0;
	return;
}

void Orbit::set_fppb_vrms(int step, float fv) {
	fppb_vrms[step] = fv;
	if (fv > 0)
		fppb[step] = 1;
	if (fv <= 0)
		fppb[step] = 0;
	return;
}

void Orbit::set_fppb_eq_snp(int step, bool x) {
	fppb_eq_snp[step] = x;
	if (fppb_mass[step] != snp_mass[step])
		cout << "WARNING: FPPB_EQ_SNP flag set with FPPB mass != SNP mass."
				<< endl;
	return;
}

void Orbit::set_host_id(int step, int id) {
	host_id[step] = id;
	return;
}

void Orbit::set_xp(int step, float xval) {
	xp[step] = xval;
	return;
}

void Orbit::set_yp(int step, float yval) {
	yp[step] = yval;
	return;
}

void Orbit::set_zp(int step, float zval) {
	zp[step] = zval;
	return;
}

void Orbit::set_vxp(int step, float vxval) {
	vxp[step] = vxval;
	return;
}

void Orbit::set_vyp(int step, float vyval) {
	vyp[step] = vyval;
	return;
}

void Orbit::set_vzp(int step, float vzval) {
	vzp[step] = vzval;
	return;
}

void Orbit::set_xf(int step, float xval) {
	xf[step] = xval;
	return;
}

void Orbit::set_yf(int step, float yval) {
	yf[step] = yval;
	return;
}

void Orbit::set_zf(int step, float zval) {
	zf[step] = zval;
	return;
}

void Orbit::set_vxf(int step, float vxval) {
	vxf[step] = vxval;
	return;
}

void Orbit::set_vyf(int step, float vyval) {
	vyf[step] = vyval;
	return;
}

void Orbit::set_vzf(int step, float vzval) {
	vzf[step] = vzval;
	return;
}

//other functions

void Orbit::display() {
	cout << "#Orbit Header: size_fppb, branch_length" << endl;
	cout
			<< "#Snap data line 1: ID, Scale, Halo Mass, Parent Mass, Subsub, Super^n Parent (SNP) Mass, Final Parent Progenitor Branch (FPPB) Mass, FPPB=SNP?, SNP ID"
			<< endl;
	cout << "#Snap data line 2: Relative X (to parent), Y, Z, VX, VY, VZ"
			<< endl;
	cout << "#Snap data line 3: Relative X (to FPPB), Y, Z, VX, VY, VZ" << endl;

	cout << nsteps << endl;
	for (int loop = 0; loop < nsteps; loop++) {
		cout << id[loop] << " " << scale[loop] << " " << mass[loop] << " "
				<< parent_mass[loop] << " " << subsub[loop] << " "
				<< snp_mass[loop] << " " << fppb_mass[loop] << " "
				<< fppb_eq_snp[loop] << " " << endl;
		if (parent[loop] == 1) {
			cout << xp[loop] << " " << yp[loop] << " " << zp[loop] << " "
					<< vxp[loop] << " " << vyp[loop] << " " << vzp[loop] << " "
					<< endl;
		} else {
			cout << "VOID" << endl;
		}
		if (fppb[loop] == 1) {
			cout << xf[loop] << " " << yf[loop] << " " << zf[loop] << " "
					<< vxf[loop] << " " << vyf[loop] << " " << vzf[loop] << " "
					<< endl;
		} else {
			cout << "VOID" << endl;
		}
	}
	return;
}

float Orbit::rp(int step) {
	return sqrt(xp[step] * xp[step] + yp[step] * yp[step] + zp[step] * zp[step]);
}

float Orbit::vp(int step) {
	return sqrt(
			vxp[step] * vxp[step] + vyp[step] * vyp[step]
					+ vzp[step] * vzp[step]);
}

float Orbit::rf(int step) {
	return sqrt(xf[step] * xf[step] + yf[step] * yf[step] + zf[step] * zf[step]);
}

float Orbit::vf(int step) {
	return sqrt(
			vxf[step] * vxf[step] + vyf[step] * vyf[step]
					+ vzf[step] * vzf[step]);
}

float Orbit::vrp(int step) {
	return (vxp[step] * xp[step] + vyp[step] * yp[step] + vzp[step] * zp[step])
			/ sqrt(
					xp[step] * xp[step] + yp[step] * yp[step]
							+ zp[step] * zp[step]);
}

float Orbit::vrf(int step) {
	return (vxf[step] * xf[step] + vyf[step] * yf[step] + vzf[step] * zf[step])
			/ sqrt(
					xf[step] * xf[step] + yf[step] * yf[step]
							+ zf[step] * zf[step]);
}

float Orbit::vtp(int step) {
	return sqrt(
			(vxp[step] * vxp[step] + vyp[step] * vyp[step]
					+ vzp[step] * vzp[step]) - vrp(step) * vrp(step));
}

float Orbit::vtf(int step) {
	return sqrt(
			(vxf[step] * vxf[step] + vyf[step] * vyf[step]
					+ vzf[step] * vzf[step]) - vrf(step) * vrf(step));
}

float Orbit::rpproj(int step, int proj) {
	if (proj == 0) {
		return sqrt(xp[step] * xp[step] + yp[step] * yp[step]);
	}
	if (proj == 1) {
		return sqrt(yp[step] * yp[step] + zp[step] * zp[step]);
	}
	if (proj == 2) {
		return sqrt(zp[step] * zp[step] + xp[step] * xp[step]);
	} else {
		cout << "Invalid projection (radius)." << endl;
		exit(0);
	}
}

float Orbit::rfproj(int step, int proj) {
	if (proj == 0) {
		return sqrt(xf[step] * xf[step] + yf[step] * yf[step]);
	}
	if (proj == 1) {
		return sqrt(yf[step] * yf[step] + zf[step] * zf[step]);
	}
	if (proj == 2) {
		return sqrt(zf[step] * zf[step] + xf[step] * xf[step]);
	} else {
		cout << "Invalid projection (radius)." << endl;
		exit(0);
	}
}

float Orbit::vpproj(int step, int proj) {
	if (proj == 0) {
		return vzp[step]
				+ brgastro::H_0 * zp[step] * snp_radius[step] / (10.0 * snp_vrms[step]);
	}
	if (proj == 1) {
		return vxp[step]
				+ brgastro::H_0 * xp[step] * snp_radius[step] / (10.0 * snp_vrms[step]);
	}
	if (proj == 2) {
		return vyp[step]
				+ brgastro::H_0 * yp[step] * snp_radius[step] / (10.0 * snp_vrms[step]);
	} else {
		cout << "Invalid projection (velocity)." << endl;
		exit(0);
	}
}

float Orbit::vfproj(int step, int proj) {
	if (proj == 0) {
		return vzf[step]
				+ brgastro::H_0 * zf[step] * fppb_radius[step] / (10.0 * fppb_vrms[step]);
	}
	if (proj == 1) {
		return vxf[step]
				+ brgastro::H_0 * xf[step] * fppb_radius[step] / (10.0 * fppb_vrms[step]);
	}
	if (proj == 2) {
		return vyf[step]
				+ brgastro::H_0 * yf[step] * fppb_radius[step] / (10.0 * fppb_vrms[step]);
	} else {
		cout << "Invalid projection (velocity)." << endl;
		exit(0);
	}
}

#if (1) //brg

brgastro::stripping_orbit *Orbit::stripping_orbit_spline() {
	if (!orbit_spline_loaded) throw std::runtime_error("Must load orbit spline before accessing it.");
	return orbit_spline_ptr;
}

const int Orbit::load_orbit_spline(const int init_step, brgastro::stripping_orbit *new_ptr,
		const bool silent)
{
	if (orbit_spline_loaded) {
		return 0; // Already loaded
	}

	int rval = 0;

#pragma omp critical
	{
		double min_good_d = 1e-5;
		double max_good_d = 1e5;
		double min_good_v = 1e-5;
		double max_good_v = 1e5;

		bool first_good_snap = true;

		// Not already loaded, so we'll load it

		try {

			orbit_spline_ptr = new_ptr;
			orbit_spline_ptr->clear();

			double init_mass = 1;
			double last_t = (-DBL_MAX);
			int last_host_id = -1;

			for (int loop = init_step; loop < nsteps; loop++) {
				bool good_snap = true;
				// halo mass
				if (mass.at(loop) <= 0) {
					good_snap = false;
				}
				//snp mass
				if (snp_mass.at(loop) <= 0) {
					good_snap = false;
				}
				//snp radius
				if (snp_radius.at(loop) == -1) {
					good_snap = false;
				}
				//snp vrms
				if (snp_vrms.at(loop) == -1) {
					good_snap = false;
				}
				if (xp.at(loop) == -1) {
					good_snap = false;
				}

				if (good_snap) {
					// It seems good so far, so now check if positions are good or glitched
					if ((std::fabs(xp.at(loop)) < min_good_d)
							|| (std::fabs(yp.at(loop)) < min_good_d)
							|| (std::fabs(zp.at(loop)) < min_good_d)
							|| (std::fabs(xp.at(loop)) > max_good_d)
							|| (std::fabs(yp.at(loop)) > max_good_d)
							|| (std::fabs(zp.at(loop)) > max_good_d)) {
						good_snap = false;
					}
				}

				double current_t = brgastro::tfa(scale.at(loop));
				if(current_t==last_t) good_snap = false;

				if (good_snap) {
					double current_snp_mass = snp_mass.at(loop)
							* unitconv::Msuntokg;
					int current_host_id = host_id.at(loop);
					double dfactor = snp_radius.at(loop) * scale.at(loop)
							* unitconv::kpctom;
					double vfactor = snp_vrms.at(loop) * unitconv::kmpstomps;

					if (first_good_snap) {
						// This is the beginning of the orbit, so load the initial satellite and parent masses
						init_mass = mass.at(loop) * unitconv::Msuntokg;
						orbit_spline_ptr->set_tNFW_init_satellite(init_mass,
								brgastro::zfa(scale.at(loop)));
						orbit_spline_ptr->set_tNFW_init_host(current_snp_mass,
								brgastro::zfa(scale.at(loop)));
						first_good_snap = false;
					} else {
						// Check for a discontinuity
						if (current_host_id != last_host_id) {
							orbit_spline_ptr->add_discontinuity_time(
									(current_t + last_t) / 2);
						}
					}
					// Check if the velocities are good or glitched
					if ((std::fabs(vxp.at(loop)) < min_good_v)
							|| (std::fabs(vyp.at(loop)) < min_good_v)
							|| (std::fabs(vzp.at(loop)) < min_good_v)
							|| (std::fabs(vxp.at(loop)) > max_good_v)
							|| (std::fabs(vyp.at(loop)) > max_good_v)
							|| (std::fabs(vzp.at(loop)) > max_good_v)) {
						orbit_spline_ptr->add_point(xp.at(loop) * dfactor,
								yp.at(loop) * dfactor, zp.at(loop) * dfactor,
								current_t,
								mass.at(loop) * unitconv::Msuntokg
										/ brgastro::safe_d(init_mass));
					} else {
						orbit_spline_ptr->add_point(xp.at(loop) * dfactor,
								yp.at(loop) * dfactor, zp.at(loop) * dfactor,
								vxp.at(loop) * vfactor, vyp.at(loop) * vfactor,
								vzp.at(loop) * vfactor, current_t,
								mass.at(loop) * unitconv::Msuntokg
										/ brgastro::safe_d(init_mass));
					}
					unsigned int num_host_parameters = 4;
					std::vector<BRG_UNITS> host_parameters(num_host_parameters, 0);
					try {
						host_parameters.at(0) = current_snp_mass;
						host_parameters.at(1) = brgastro::zfa(scale.at(loop));
						host_parameters.at(2) = 0; // Use default c
						host_parameters.at(3) = 0; // Use default tau

						orbit_spline_ptr->add_host_parameter_point(host_parameters, current_t);
					} catch (...) {
						if(!silent) std::cerr
								<< "WARNING: Could not set host parameters in orbit.\n";
					}

					last_host_id = current_host_id;
					last_t = current_t;
				}
			}

		} catch (const std::exception &) {
			std::cerr << "ERROR: Could not read in orbit spline.\n";
			orbit_spline_loaded = false;
			rval = UNSPECIFIED_ERROR;
		}

		if (first_good_snap) {
			// If this is the case, we were unable to find any good snapshots. Don't load any orbit for it
			orbit_spline_loaded = false;
			rval = UNSPECIFIED_ERROR;
		}

		orbit_spline_loaded = true;

	}

	return rval;

}

void Orbit::unload_orbit_spline() {
	orbit_spline_loaded = false;
}
#endif //brg

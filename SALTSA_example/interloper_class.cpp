#include "interloper_class.h"

using namespace std;

Interloper::Interloper() {
	ID = -1;
	mass = -1;
	host_mass = -1;
	host_rvir = -1;
	host_vrms = -1;
	R = 0;
	V = 0;
	return;
}

Interloper::Interloper(ifstream& data) {
	// Backup initialisation in case we can't read in data
	ID = -1;
	mass = -1;
	host_mass = -1;
	host_rvir = -1;
	host_vrms = -1;
	R = 0;
	V = 0;

	string junk;
	data >> ID;
	data >> mass;
	data >> host_mass;
	data >> host_rvir;
	data >> host_vrms;
	data >> R;
	data >> V;
	getline(data, junk);
	return;
}

Interloper::Interloper(const Interloper& other) {
	ID = other.ID;
	mass = other.mass;
	host_mass = other.host_mass;
	host_rvir = other.host_rvir;
	host_vrms = other.host_vrms;
	R = other.R;
	V = other.V;
	return;
}

void Interloper::swap(Interloper& other) {
	std::swap(ID, other.ID);
	std::swap(mass, other.mass);
	std::swap(host_mass, other.host_mass);
	std::swap(host_rvir, other.host_rvir);
	std::swap(host_vrms, other.host_vrms);
	std::swap(R, other.R);
	std::swap(V, other.V);
	return;
}

#ifndef SWIG
Interloper& Interloper::operator=(Interloper other) {
	Interloper::swap(other);
	return *this;
}
#endif

Interloper::~Interloper() {
	return;
}

long int Interloper::get_ID() {
	return ID;
}

float Interloper::get_mass() {
	return mass;
}

float Interloper::get_host_mass() {
	return host_mass;
}

float Interloper::get_host_rvir() {
	return host_rvir;
}

float Interloper::get_host_vrms() {
	return host_vrms;
}

float Interloper::get_R() {
	return R;
}

float Interloper::get_V() {
	return V;
}

void Interloper::set_ID(long int id) {
	ID = id;
	return;
}

void Interloper::set_mass(float m) {
	mass = m;
	return;
}

void Interloper::set_host_mass(float hm) {
	host_mass = hm;
	return;
}

void Interloper::set_host_rvir(float hr) {
	host_rvir = hr;
	return;
}

void Interloper::set_host_vrms(float hv) {
	host_vrms = hv;
	return;
}

void Interloper::set_R(float r) {
	R = r;
	return;
}

void Interloper::set_V(float v) {
	V = v;
	return;
}

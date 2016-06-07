#include "read_orbits.h"
#include <limits>

using namespace std;

#ifdef ORB_LIST
void read_orbit_data(const char* datafile, list <Orbit>& orbits,
		unsigned int N_start, unsigned int N_to_read)
#else
void read_orbit_data(const char* datafile, vector<Orbit>& orbits,
		unsigned int N_start, unsigned int N_to_read)
#endif
		{
	orbits.clear();
	ifstream data;
	data.open(datafile, ios::in);

	string junk;

	if(N_to_read==0) N_to_read = std::numeric_limits<unsigned int>::max()-N_start;

	if (data.is_open()) {
		unsigned int N = 0;
		while (!data.eof()) {
			if (data.peek() == '#') {
				getline(data, junk);
				continue;
			}
			Orbit tmp = Orbit(data);
			if( (N >= N_start) && (N < N_start+N_to_read) )
				orbits.push_back(tmp);
			++N;
			//getline(data, junk);
			if (data.peek() < 0)
				break;
		}
	} else {
		cout << "Failed to open datafile." << endl;
		exit(0);
	}
}

void read_interloper_data(const char* datafile,
		vector<Interloper>& interlopers) {
	interlopers.clear();
	ifstream data;

	data.open(datafile, ios::in);

	string junk;

	if (data.is_open()) {
		while (!data.eof()) {
			if (data.peek() == '#') {
				getline(data, junk);
				continue;
			}
			interlopers.push_back(Interloper(data));
			if (data.peek() < 0)
				break;
		}
	} else {
		cout << "Failed to open datafile." << endl;
		exit(0);
	}
	return;
}

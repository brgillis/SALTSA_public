#! /bin/bash
echo "Attempting to compile SALTSA_example..."
make
make clean
echo "Attempting to run SALTSA_example..."
./SALTSA_example
echo "Attempting to generate orbit plot (requires perl and stilts)..."
./make_SALTSA_example_plot.pl
echo "Attempting to display orbit plot using gpicview..."
gpicview SALTSA_example_orbits.png

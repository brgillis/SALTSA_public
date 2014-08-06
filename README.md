### What is this repository for? ###

* Semi-AnaLytic Tidal Stripping Algorithm (SALTSA)
* Version 1.0

### How do I get set up? ###

* Dependencies
    * A half-decent C++ compiler (doesn't even need to handle C++0x/C++11).

* How to get set up
    * Download the repository, which includes an example script

* How to run tests
    * On a Linux-like system which can run bash scripts, run `./EXECUTEME.sh` to do a quick compile and test
        * This script will try to use `make` to compile. The current set-up assumes use of g++ for the compiler. You'll have to tweak things for another compiler.
        * The script also tries to display a plot of the results. This part requires perl and stilts to be installed. If it doesn't work, don't worry, that's just one way to display the results.
        * If it does work, the displayed image should match the below:

        ![What your test output from SALTSA_example with the default settings should look like](https://bitbucket.org/brgillis/saltsa_public/raw/master/SALTSA_example/SALTSA_example_plot.png)

* How to use it
    * Add all *.cpp files to any C++ project except SALTSA_example.cpp. In any file that wants to access SALTSA, put the line `#include "SALTSA.h"` at the top of it, and then follow the example script in how it works.
    * You can optionally remove the "SALTSA_example.cpp" file and compile SALTSA as a library, and then link to it from your project.

### Who do I talk to? ###

* Contact Dr Bryan Gillis if you find any bugs in the code or you wish to discuss extensions to it.
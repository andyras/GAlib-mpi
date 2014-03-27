#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <numeric>
#include <complex>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_diag.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>
#include <map>
#include <fftw3.h>
#include <omp.h>

#include "libdynamix_input_parser.hpp"
#include "libdynamix_outputs.hpp"
#include "output.hpp"
#include "numerical.hpp"
#include "params.hpp"
#include "userdata.hpp"
#include "rhs.hpp"
#include "plots.hpp"
#include "constants.hpp"
#include "conversions.hpp"
#include "analytic.hpp"

/* main function for dynamix */
int dynamixMain (int argc, char * argv[]);

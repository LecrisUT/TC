#ifndef TC_HEADER_HPP
#define TC_HEADER_HPP

#define NDEBUG

// all quantities are in atomic unit in TC++

#include <cmath>
#include <cstdlib>

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <complex>
#include <vector>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "fftw3.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/property_tree/xml_parser.hpp>

using Complex = std::complex<double>;
constexpr double PI = 3.14159265358979323846; // std::numbers::pi in C++20
constexpr double FourPI = 4 * PI;
// constexpr Complex I {0.0,1.0}; // constexpr Complex does not work for some
// (old) compilers...
const Complex I{0.0, 1.0};
constexpr double Ht_in_eV = 27.21138505;

// TC++ header files
#include "error_messages.hpp"
#include "file_names.hpp"
#include "method.hpp"
#include "my_clock.hpp"
#include "spin.hpp"

#include "crystal_structure.hpp"
#include "kpoints.hpp"
#include "symmetry.hpp"

#include "bloch_states.hpp"
#include "plane_wave_basis.hpp"
#include "total_energy.hpp"

#include "parallelization.hpp"

#include "calc_hamiltonian.hpp"
#include "diagonalization.hpp"
#include "jastrow.hpp"
#include "potentials.hpp"

#include "io_tc_files.hpp"
#include "read_qe.hpp"
#include "read_upf.hpp"

#endif // TC_HEADER_HPP

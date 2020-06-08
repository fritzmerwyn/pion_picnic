#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include "help_numerics.h"
// #include "progressbar.hpp"

std::complex<double> funccomplex(std::complex<double> z);

int precalculation(double m_c, double* renorm_constants, double* a_vals, double* b_vals, double*absciss_x, double*absciss_ang, double*weights_w, double*weights_ang, double eta);

#ifndef help_H    // So you dont declare header/function twice
#define help_H

// Gauss Legendre Code from Numerical Recipes in C
// and more numerical functions for the integration.
#include "definitions.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#define EPS 3.0e-11 // EPS is the relative precision.

double** gauleg(double a,double b,int e);
double angkern(double x);
double angkern2(double x);
double qgaus1(double (*func)(double), double* x, double* w);
int locate(double* xx, double x, int length);
void read_in_data(const std::string filename, double* q_vec, double* z_vec, std::complex<double>*** y_vals);
double* get_qz(std::complex<double> x0, double m_pion, double routing_plus);
std::complex<double> bilinearinterpol(double q, double z, double* x_corners, std::complex<double>* y_corners);

#endif

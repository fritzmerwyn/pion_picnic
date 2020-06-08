#include <math.h>
#include <iostream>
#include "help_numerics.h"
#include <complex>
#include <fstream>
#define EPS 3.0e-11 // EPS is the relative precision.

double qgaus1(double (*func)(double), double* x, double* w){
  int j;
  double s;
  s=0.0;
  for(j=0;j<=absciss_points;j++){
    s += w[j]*((*func)(x[j]));
  }
  return s;
}

double angkern(double x){
  return sqrt(1-x*x);
}

double angkern2(double x){
  return sqrt(1-x*x)*x;
}

double** gauleg(double x1, double x2, int n)
// Given the lower and upper limits of integration x1 and x2, and given n,
// this routine returns arrays x[1..n] and w[1..n] of length n, containing
// the abscissas and weights of the Gauss- Legendre n-point quadrature formula.
{
int m,j,i;
double z1,z,xm,xl,pp,p3,p2,p1;
double* x = nullptr;
double* w = nullptr;

x = new double[n+1]; // WARNING the last entry (n+1) is still there but is the same as n. Make sure the other programs do nut run until the last entry <=n but only <n!
w = new double[n+1]; // WARNING the last entry (n+1) is still there but is the same as n. Make sure the other programs do nut run until the last entry <=n but only <n!

m=(n+1)/2;
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);
for(i=1;i<=m;i++) {
  z=cos(M_PI*(i-0.25)/(n+0.5));
  // Starting with the above approximation to the ith root, we enter the main loop of refinement by Newton’s method.
  do {
    p1=1.0;
    p2=0.0;
    for (j=1;j<=n;j++) {
      p3=p2;
      p2=p1;
      p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
    };
    // p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2, the polynomial of one lower order.
    pp=n*(z*p1-p2)/(z*z-1.0);
    z1=z;
    z=z1-p1/pp;
  } while (fabs(z-z1) > EPS);
  x[i]=xm-xl*z;
  x[n+1-i]=xm+xl*z;
  w[i]=2.0*xl/((1.0-z*z)*pp*pp);
  w[n+1-i]=w[i];
// Scale the root to the desired interval, and put in its symmetric counterpart. Compute the weight
// and its symmetric counterpart.
  }

  for(int k = 0;k<n; k++){
    x[k] = x[k+1]; // WARNING the last entry (n+1) is still there but is the same as n. Make sure the other programs do nut run until the last entry <=n but only <n!
    w[k] = w[k+1]; // WARNING the last entry (n+1) is still there but is the same as n. Make sure the other programs do nut run until the last entry <=n but only <n!
  }



  // x.pop();
  // w.pop();
  // Initialises a pointer to a 2d array (abcissas_gl2d) which contains abscissas x[j] and weights
  // w[j] in [0] and [1] respectively.
  // This array is then returned.
  double** abcissas_gl2d = nullptr;
  abcissas_gl2d = new double*[2];
  abcissas_gl2d[0] = x;
  abcissas_gl2d[1] = w;

  return abcissas_gl2d;
}

int locate(double* xx, double x, int length){

    unsigned long ju,jm,jl;
    int ascnd;

    if(x == xx[0]) return 0;
    else if(x == xx[length-1]) return length-2;

    jl=0;
    ju=length;
    ascnd=(xx[length-1] >= xx[0]);
    while(ju-jl > 1){
      // jm=(ju+jl) >> 1; // original code... shifts the bit value of jm one bit to right eg from 6=00110 to 3=00011
      jm=(ju+jl)/2; // original code... shifts the bit value of jm one bit to right eg from 6=00110 to 3=00011

      if (x >= xx[jm] == ascnd){
        jl=jm;
      }
      else{
        ju=jm;
      }
    }

  return jl;
}

void read_in_data(const std::string filename, double* q_vec, double* z_vec, std::complex<double>*** y_vals){
  std::complex<double> Imag = {0.0,1.0};
  double x0,x1,x2,x3,x4,x5;
  int i,j,k;
  j=0;
  k=0;
  std::ifstream inputfile(filename);
  inputfile.ignore(50, '\n');

  if (!inputfile) {
    std::cerr << "Unable to open specified file... Quark already iterated?";
    exit(1);   // call system to stop
  }

  while (inputfile >> x0 >> x1 >> x2 >> x3 >> x4 >> x5) {

      if(k == 0){
        q_vec[j] = x0;
      }
      if(j == 0){
        z_vec[k] = x1;
      }

      y_vals[0][j][k] = x2 + Imag*x3; // A+
      y_vals[1][j][k] = x4 + Imag*x5; // B+

      k++;
      if(k == ang_absciss_points){
        k=0;
        j++;
      }
  }
  inputfile.close();

}

double* get_qz(std::complex<double> x0, double m_pion, double routing_plus){
  double q = std::sqrt(x0.real() + routing_plus * routing_plus * m_pion * m_pion);
  double z = (x0.imag())/(2.0*routing_plus * m_pion * q);

  double* got_qz = nullptr;
  got_qz = new double[2];
  got_qz[0] = q;
  got_qz[1] = z;

  return got_qz;
}

std::complex<double> bilinearinterpol(double q, double z, double* x_corners, std::complex<double>* y_corners){
  double t,u;

  // double* qz = nullptr;
  //
  // qz = get_qz(x0,m_pion,routing_plus);

  t = (q - x_corners[0])/(x_corners[1]-x_corners[0]);
  u = (z - x_corners[2])/(x_corners[3]-x_corners[2]);

  return (1-t)*(1-u)*y_corners[0] + t*(1-u)*y_corners[1] + t*u*y_corners[2] + (1-t)*u*y_corners[3];
}

#include "lookuptable.h"
#include "Dyson_test.h"
#include "progressbar.hpp"


std::complex<double> funccomplex(std::complex<double> z){
  return z*z*z - 2.0*z*z + 3.0*z - 4.0;
}

int precalculation(double m_c, double* renorm_constants, double* a_vals, double* b_vals, double*absciss_x, double*absciss_ang, double*weights_w, double*weights_ang, double eta){

  // std::complex<double>** temp_vals = nullptr;
  // temp_vals = new std::complex<double>*[3];
  // temp_vals[0] = nullptr;
  // temp_vals[1] = nullptr;
  // temp_vals[2] = nullptr;
  // temp_vals[0] = new std::complex<double>[absciss_points];
  // temp_vals[1] = new std::complex<double>[absciss_points];
  // temp_vals[2] = new std::complex<double>[absciss_points];

  double routing_plus, routing_minus, q, z, psi, m_pion;
  std::complex<double> x_c, cmplxval,x_c_sqrt;

  m_pion = PI_MAX;
  std::complex<double> Imag = {0.0,1.0};
  routing_plus = 0.5;
  routing_minus = routing_plus - 1.0;

  // std::ofstream  fileouta;
  // fileouta.open("Data/LOOKUPTABLE1_fine.dat");
  // fileouta<<"#q"<< "\t" << "z"<<"\t"<<"cmplxval.real()"<<"\t"<<"cmplxval.imag()"<<std::endl;
  //
  // for(int q_idx = 0; q_idx < absciss_points; q_idx++){
  //   q = exp(0.5*absciss_x[q_idx]);
  //   for(int psi_idx=0; psi_idx < ang_absciss_points; psi_idx++){
  //       psi = absciss_ang[psi_idx];
  //       z = std::cos(psi);
  //       // x_c = std::pow(q,2.0) + routing_minus*q*Imag*m_pion*z + routing_plus*q*Imag*m_pion*z - routing_plus*routing_minus*m_pion*m_pion;
  //       x_c = std::pow(q,2.0) + 2.0*routing_plus*q*Imag*m_pion*z - std::pow(routing_plus,2.0)*m_pion*m_pion;
  //       x_c_sqrt = std::sqrt(x_c);
  //       // cmplxval = weights_w[q_idx]*weights_ang[psi_idx]*funccomplex(x_c);
  //       cmplxval = funccomplex(x_c_sqrt);
  //
  //       fileouta<< q << "\t" << z <<"\t"<< 0.0 <<"\t" << 0.0 << "\t" << cmplxval.real()<<" "<<cmplxval.imag()<<std::endl;
  //     }
  //   }
  // fileouta.close();

  std::ofstream  fileoutb;
  fileoutb.open("Data/DressingFunctionsABCmplxInterpol_even_rougherPI_MAX.dat");
  fileoutb<<"#q"<< "\t" << "z" <<"\t" <<"A.real()"<< "\t" << "A.imag()"<<"\t"<<"B.real()"<<"\t"<<"B.imag()"<<std::endl;

  ProgressBar qi(absciss_points, "Iterating Quark");

  for(int q_idx = 0; q_idx < absciss_points; q_idx++){
    ++qi;
    q = exp(0.5*absciss_x[q_idx]);
    for(int psi_idx=0; psi_idx < ang_absciss_points; psi_idx++){
        psi = absciss_ang[psi_idx];
        z = std::cos(psi);
        // x_c = std::pow(q,2.0) + routing_minus*q*Imag*m_pion*z + routing_plus*q*Imag*m_pion*z - routing_plus*routing_minus*m_pion*m_pion;
        x_c = std::pow(q,2.0) + 2.0*routing_plus*q*Imag*m_pion*z - std::pow(routing_plus,2.0)*m_pion*m_pion;
        x_c_sqrt = std::sqrt(x_c);
        // cmplxval = weights_w[q_idx]*weights_ang[psi_idx]*funccomplex(x_c);
        std::complex<double>* interpolvalsAB = interpolation_cmplx(x_c_sqrt, m_c, renorm_constants, a_vals, b_vals, absciss_x, absciss_ang, weights_w, weights_ang, eta);


        // temp_vals[0][psi_idx] = x_c;
        // temp_vals[1][psi_idx] = cmplxval.real();
        // temp_vals[2][psi_idx] = cmplxval.imag();

        fileoutb<< q << " " << z << " " << interpolvalsAB[0].real() << " " << interpolvalsAB[0].imag() << " " << interpolvalsAB[1].real() << " " << interpolvalsAB[1].imag()<<std::endl;
      }
    }
  fileoutb.close();
  std::cout<<"Precalc Done!"<<std::endl;
  //
  return 0;

}

#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
//#include "kernel.h"
#include "help_numerics.h"
#include "Dyson_test.h"
#include "progressbar.hpp"

double gamma_fun(double color = N_C, double flavor = N_F){
  return 12.0/(11.0 * color - 2.0 * flavor);
}

double running_coupling_MarisTandy(double k_squared, double eta){
  double infrared, ultraviolet;
  infrared = M_PI * ((pow(eta,7.0))/(pow(Lambda_0,4.0))) * (k_squared * k_squared) * exp( (-(eta*eta)) * (k_squared/(Lambda_0*Lambda_0)) );
  ultraviolet = (2.0 * M_PI * gamma_fun() * (1.0 - exp(-(k_squared)/(Lambda_t*Lambda_t)) ) ) / (log( M_E*M_E - 1.0 + pow(1.0 + (k_squared/(Lambda_QCD*Lambda_QCD) ),2.0) ));
  return  infrared + ultraviolet;
}

std::complex<double> running_coupling_MarisTandy_cmplx(std::complex<double> k_squared, double eta){
  std::complex<double> infrared, ultraviolet;
  infrared = M_PI * ((std::pow(eta,7.0))/(std::pow(Lambda_0,4.0))) * (k_squared * k_squared) * std::exp( (-(eta*eta)) * (k_squared/(Lambda_0*Lambda_0)) );
  ultraviolet = (2.0 * M_PI * gamma_fun() * (1.0 - std::exp(-(k_squared)/(Lambda_t*Lambda_t)) ) ) / (std::log( M_E*M_E - 1.0 + std::pow(1.0 + (k_squared/(Lambda_QCD*Lambda_QCD) ),2.0) ));
  return  infrared + ultraviolet;
}

double*** initialize_matrix(double epsilon, double m_c, double* absciss_x, double* absciss_ang, double* weights_ang, double g_squared, double eta, double mu){

  double q, c2_a, c2_b, p, z, yota, k_squared, s0_a, s0_b;

  c2_a = 4.0/(2.0*3.0*pow(M_PI,2.0));
  c2_b = 4.0/(2.0*pow(M_PI,2.0));


  double*** temp_matrix = nullptr;
  temp_matrix = new double**[2];
  temp_matrix[0] = nullptr;
  temp_matrix[1] = nullptr;
  temp_matrix[0] = new double*[absciss_points];
  temp_matrix[1] = new double*[absciss_points];

#pragma omp parallel for default(none) shared(temp_matrix)
  for(int i=0;i<absciss_points;i++){
    temp_matrix[0][i] = nullptr;
    temp_matrix[1][i] = nullptr;
    temp_matrix[0][i] = new double[absciss_points + 1]; // + 1, so that in the last index, the value for p=mu is saved.
    temp_matrix[1][i] = new double[absciss_points + 1];
  }
  ProgressBar initmatrix(absciss_points,"Initializing Angular Matrix: ");
  for(int q_idx = 0; q_idx < absciss_points; q_idx++, ++initmatrix){

    q = exp(0.5*absciss_x[q_idx]);

    for(int p_idx = 0; p_idx < absciss_points + 1; p_idx++){

      if(p_idx == absciss_points){
        p = mu;
      }

      else{
      p = exp(0.5*absciss_x[p_idx]);
      }

      s0_a=0.0;
      s0_b=0.0;

#pragma omp parallel for private(z, yota, k_squared) default(none) shared(p, q, weights_ang, absciss_ang, absciss_x, c2_a, c2_b, eta, mu) reduction(+:s0_a, s0_b)
            for(int ang_idx=0;ang_idx<ang_absciss_points;ang_idx++){ //START AT J=1 because 0th abscissa is 0. And 0 no good.(look at bottom of main)
              yota = absciss_ang[ang_idx];
              z = cos(yota);
              k_squared = p*p + q*q - 2.0*p*q*z;
              if(p==0.0){
              }
              else{
              s0_a += (c2_a * q*q*q*q)* (1.0/(p*p)) *
                      weights_ang[ang_idx] * sin(yota)*sin(yota) * (p*q*z + (2.0/(k_squared)) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z))   *
                      (running_coupling_MarisTandy(k_squared, eta) / (k_squared));

              s0_b += (c2_b * q*q*q*q) *
                      weights_ang[ang_idx] * (sin(yota)*sin(yota) * running_coupling_MarisTandy(k_squared,eta)/(k_squared));
              }
            }

      std::swap(temp_matrix[0][q_idx][p_idx],s0_a);
      std::swap(temp_matrix[1][q_idx][p_idx],s0_b);
      // temp_matrix[0][q_idx][p_idx] = s0_a;
      // temp_matrix[1][q_idx][p_idx] = s0_b;

    }
  }

  return temp_matrix;

  // // ##### FREE matrix pointer
  // for (int i = 0; i<absciss_points ; i++){
  //   free(temp_matrix[0][i]);
  //   free(temp_matrix[1][i]);
  // }
  // free(temp_matrix[0]);
  // free(temp_matrix[1]);
  // free(temp_matrix);

}

double** initialize_dressing_functionAB(double a0, double b0){
  double** dress2d = nullptr;
  dress2d = new double*[3];
  dress2d[0] = nullptr;
  dress2d[1] = nullptr;
  dress2d[2] = nullptr;
  dress2d[0] = new double[absciss_points];
  dress2d[1] = new double[absciss_points];
  dress2d[2] = new double[2]; // 2 is the number of renormalization constants used in the iteration-function

  ProgressBar idf(absciss_points, "Initializing Dressing Functions");

  for(int i=0;i<absciss_points;i++,++idf){
    dress2d[0][i] = a0;
    dress2d[1][i] = b0;
  }
    dress2d[2][0] = 1.0;
    dress2d[2][1] = 1.0;

  std::cout<<std::endl;
  std::cout<<"Dressing Functions initialized as A = "<< a0 <<" and B = "<< b0 << " z2 as: " << dress2d[2][0]<< " zm as: "<< dress2d[2][1]<< std::endl;
  return dress2d;

  // ##### FREE Vector pointer #####
  // for (int i = 0; i<3 ; i++){
  //   free(dress2d[i]);
  // }
  // free(dress2d);

}

std::complex <double>***** initialize_theta_matrix(double* renorm_constants,double* BSE_absciss_x, double* BSE_absciss_ang, double* BSE_weights_w,
                                                  double* BSE_weights_ang, double eta, double alpha){

  std::complex<double>***** temp_mother_matrix_1 = nullptr;
  temp_mother_matrix_1 = new std::complex<double>****[2];
  // std::complex<double>*** temp_mother_matrix_2 = nullptr;
  // temp_mother_matrix_2 = new std::complex<double>**[BSE_ang_absciss_points];
  // std::complex<double>*** temp_mother_matrix_3 = nullptr;
  // temp_mother_matrix_3 = new std::complex<double>**[BSE_ang_absciss_points];
  // std::complex<double>*** temp_mother_matrix_4 = nullptr;
  // temp_mother_matrix_4 = new std::complex<double>**[BSE_ang_absciss_points];

  for(int i=0;i<2;i++){
    temp_mother_matrix_1[i] = nullptr;
    temp_mother_matrix_1[i] = new std::complex<double>***[BSE_ang_absciss_points];
  }

  for(int i=0;i<2;i++){
    for(int j=0;j<BSE_ang_absciss_points;j++){
      temp_mother_matrix_1[i][j] = nullptr;
      temp_mother_matrix_1[i][j] = new std::complex<double>**[BSE_ang_absciss_points];
    }
  }

  for(int i=0;i<2;i++){
    for(int j=0;j<BSE_ang_absciss_points;j++){
      for(int k=0;k<BSE_ang_absciss_points;k++){
    temp_mother_matrix_1[i][j][k] = nullptr;
    temp_mother_matrix_1[i][j][k] = new std::complex<double>*[BSE_absciss_points];
      }
    }
    // temp_mother_matrix_2[i] = nullptr;
    // temp_mother_matrix_2[i] = new std::complex<double>*[BSE_absciss_points];
    // temp_mother_matrix_3[i] = nullptr;
    // temp_mother_matrix_3[i] = new std::complex<double>*[BSE_absciss_points];
    // temp_mother_matrix_4[i] = nullptr;
    // temp_mother_matrix_4[i] = new std::complex<double>*[BSE_absciss_points];
  }

  for(int i=0;i<2;i++){
    for(int j=0;j<BSE_ang_absciss_points;j++){
      for(int k=0;k<BSE_ang_absciss_points;k++){
        for(int l=0;l<BSE_absciss_points;l++){

    temp_mother_matrix_1[i][j][k][l] = nullptr;
    temp_mother_matrix_1[i][j][k][l] = new std::complex<double>[BSE_absciss_points];
        }
      }
    }
  }

  double q, z, psi, theta,real_thetamatrix_entry_1, imag_thetamatrix_entry_1, real_thetamatrix_entry_2, imag_thetamatrix_entry_2, real_thetamatrix_entry_3, imag_thetamatrix_entry_3,
    real_thetamatrix_entry_4, imag_thetamatrix_entry_4, realMTtemp, imagMTtemp;
  std::complex<double> matrix_entry, p, k_squared;
  std::complex<double> Imag = {0.0,1.0};



#pragma omp parallel for default(none) \
                          shared(BSE_absciss_x,BSE_absciss_ang,alpha,renorm_constants,BSE_weights_ang,BSE_weights_w,eta,temp_mother_matrix_1) \
                          private(p,q,z,psi,theta, k_squared,realMTtemp,imagMTtemp) \
                          reduction(+:real_thetamatrix_entry_1,imag_thetamatrix_entry_1)

  for(int q_idx = 0; q_idx < BSE_absciss_points; q_idx++){

    q = exp(0.5*BSE_absciss_x[q_idx]);

    for(int p_idx = 0; p_idx < BSE_absciss_points; p_idx++){

      p = exp(0.5*BSE_absciss_x[p_idx]);

      std::complex<double> matrix_entry = {0.0,0.0};
      // k_squared = p*p + q*q - 2.0*p*q*(z*cos(alpha)+sin(psi)*cos(theta)*sin(alpha));

      for(int psi_idx=0; psi_idx < BSE_ang_absciss_points; psi_idx++){

          psi = BSE_absciss_ang[psi_idx];
          z = std::cos(psi);

          imag_thetamatrix_entry_1=0.0;
          real_thetamatrix_entry_1=0.0;

          realMTtemp=0.0;
          imagMTtemp=0.0;

        for(int theta_idx=0; theta_idx < BSE_ang_absciss_points; theta_idx++){

          theta = BSE_absciss_ang[theta_idx];

          k_squared = p*p + q*q - 2.0*p*q*(z*std::cos(alpha) + std::sin(psi)*std::cos(theta)*std::sin(alpha));

          realMTtemp = ((running_coupling_MarisTandy_cmplx(k_squared,eta)/k_squared).real()) *((4.0/3.0)*4.0*M_PI*2.0*M_PI*(1.0/2.0)*std::pow((1.0/(2.0*M_PI)),4.0)*
                      std::pow(renorm_constants[0],2.0)*BSE_weights_w[q_idx]*q*q*q*q*BSE_weights_ang[psi_idx] * std::sin(psi)*std::sin(psi));
          imagMTtemp = ((running_coupling_MarisTandy_cmplx(k_squared,eta)/k_squared).imag()) *((4.0/3.0)*4.0*M_PI*2.0*M_PI*(1.0/2.0)*std::pow((1.0/(2.0*M_PI)),4.0)*
                      std::pow(renorm_constants[0],2.0)*BSE_weights_w[q_idx]*q*q*q*q*BSE_weights_ang[psi_idx] * std::sin(psi)*std::sin(psi));

          real_thetamatrix_entry_1 += (BSE_weights_ang[theta_idx]*std::sin(theta)*(realMTtemp));
          imag_thetamatrix_entry_1 += (BSE_weights_ang[theta_idx]*std::sin(theta)*(imagMTtemp));

          temp_mother_matrix_1[1][theta_idx][psi_idx][p_idx][q_idx] = {realMTtemp,imagMTtemp};

        }

        std::complex<double> thetamatrix_entry_1 = {real_thetamatrix_entry_1,imag_thetamatrix_entry_1};

        temp_mother_matrix_1[0][0][psi_idx][p_idx][q_idx] = thetamatrix_entry_1;

      }
    }
  }

  // std::cout<<"Theta-Matrix initialized"<< std::endl;

  return temp_mother_matrix_1;
}

std::complex <double>** initialize_mother_matrix(double m_pion, std::complex<double>* a_corners, std::complex<double>* b_corners, double* BSE_weights_ang,
   double* BSE_absciss_x, double* BSE_absciss_ang, std::complex<double>***** theta_matrix, double* q_vec, double* z_vec, double* x_corners, std::complex<double>*** y_corner){

  double p, q, q_q, p_p, P_P, q_p,k_k,k_q,k_p, z, psi, theta, eta, eta_minus, alpha,dressing1real,dressing1imag, dressing2real, dressing2imag, dressing3real, dressing3imag, dressing4real, dressing4imag;
  double dressing1sumreal, dressing1sumimag, dressing2sumreal, dressing2sumimag, dressing3sumreal, dressing3sumimag, dressing4sumreal, dressing4sumimag;
  std::complex<double> matrix_entry,q_plus_q_minus, q_plus_squared, q_minus_squared, k_squared, P_p, q_P,k_P,P_q;
  std::complex<double> Imag = {0.0,1.0};
  alpha = 1.7;
  int grid, alpha_idx, beta_idx;
  grid = BSE_absciss_points*BSE_ang_absciss_points*4;

  double delta_pion = PI_MAX - m_pion;

  eta = 0.5;
  eta_minus = eta - 1.0;

  std::complex<double>** mother_temp = nullptr;
  mother_temp = new std::complex<double>*[grid];

  double** kernelreal = nullptr;
  kernelreal = new double*[4];
  double** kernelimag = nullptr;
  kernelimag = new double*[4];


  for(int i=0;i<4;i++){
    kernelreal[i]=nullptr;
    kernelreal[i]= new double[4];
    kernelimag[i]=nullptr;
    kernelimag[i]= new double[4];
  }

  std::complex<double>* w = nullptr;
  w = new std::complex<double>[100];


// #pragma omp parallel for default(none) shared(mother_temp)
  for(int i=0;i<grid;i++){
    mother_temp[i] = nullptr;
    mother_temp[i] = new std::complex<double>[grid];
  }
  // std::ofstream  fileouta;
  // fileouta.open("Data/DressingFunctions_A_and_B_PLUS_complex.dat");
  // fileouta<<"qp"<<"\t"<<"A+"<<"\t"<<"B+"<<"qm"<<"\t"<<"A-"<<"\t"<<"B-"<<std::endl;
  // ProgressBar mother(BSE_absciss_points,"Initializing Angular Matrix: ");


  // #pragma omp parallel for default(none) \
  //         shared(m_pion, a_corners, b_corners, BSE_weights_ang, BSE_absciss_x, BSE_absciss_ang, theta_matrix, \
  //            q_vec, z_vec, x_corners, y_corner, alpha, Imag, delta_pion, eta, eta_minus, mother_temp, matrix_entry,w,kernelimag, kernelreal) \
  //         private(p, q, z, psi, theta, \
  //            q_plus_q_minus, q_plus_squared, q_minus_squared, q_q, p_p, P_P, P_p, P_q, k_P, k_k, k_q, k_p, q_p,\
  //            alpha_idx, beta_idx)
  for(int q_idx = 0; q_idx < BSE_absciss_points; q_idx++){
    q = exp(0.5*BSE_absciss_x[q_idx]);

    for(int p_idx = 0; p_idx < BSE_absciss_points; p_idx++){

       p = exp(0.5*BSE_absciss_x[p_idx]);

       matrix_entry = {0.0,0.0};


    for(int psi_p_idx=0; psi_p_idx < BSE_ang_absciss_points; psi_p_idx++){

      for(int psi_idx=0; psi_idx < BSE_ang_absciss_points; psi_idx++){

          psi = BSE_absciss_ang[psi_idx];
          z = std::cos(psi);


          q_plus_q_minus = std::pow(q,2.0) + eta_minus*q*Imag*m_pion*z + eta*q*Imag*m_pion*z - eta*eta_minus*m_pion*m_pion;
          q_plus_squared = std::pow(q,2.0) + 2.0*eta*q*Imag*m_pion*z - std::pow(eta,2.0)*m_pion*m_pion;
          q_minus_squared = std::pow(q,2.0) + 2.0*eta_minus*q*Imag*m_pion*z - std::pow(eta_minus,2.0)*m_pion*m_pion;
          q_q = std::pow(q,2);
          p_p = std::pow(p,2);
          P_P = (-1.0)*m_pion*m_pion;
          P_p = Imag*p*m_pion*std::cos(alpha);
          P_q = Imag*q*m_pion*z;
          k_P = P_q - P_p;

          double* g_qz = get_qz(q_plus_squared,m_pion+delta_pion,eta);

          unsigned int loc1 = locate(q_vec,g_qz[0],absciss_points);
          unsigned int loc2 = locate(z_vec,g_qz[1],ang_absciss_points);

          a_corners[0] = y_corner[0][loc1][loc2];
          a_corners[1] = y_corner[0][loc1+1][loc2];
          a_corners[2] = y_corner[0][loc1+1][loc2+1];
          a_corners[3] = y_corner[0][loc1][loc2+1];

          b_corners[0] = y_corner[1][loc1][loc2];
          b_corners[1] = y_corner[1][loc1+1][loc2];
          b_corners[2] = y_corner[1][loc1+1][loc2+1];
          b_corners[3] = y_corner[1][loc1][loc2+1];

          x_corners[0] = q_vec[loc1];
          x_corners[1] = q_vec[loc1+1];
          x_corners[2] = z_vec[loc2];
          x_corners[3] = z_vec[loc2+1];



          std::complex<double> a_plus = bilinearinterpol(g_qz[0],g_qz[1],x_corners,a_corners);
          std::complex<double> b_plus = bilinearinterpol(g_qz[0],g_qz[1],x_corners,b_corners);
          std::complex<double> a_minus = std::conj(a_plus);
          std::complex<double> b_minus = std::conj(b_plus);

          std::complex<double> sigmavp = a_plus/(q_plus_squared*a_plus*a_plus + b_plus*b_plus);
          std::complex<double> sigmavm = a_minus/(q_minus_squared*a_minus*a_minus + b_minus*b_minus);
          std::complex<double> sigmasp = b_plus/(q_plus_squared*a_plus*a_plus + b_plus*b_plus);
          std::complex<double> sigmasm = b_minus/(q_minus_squared*a_minus*a_minus + b_minus*b_minus);


          for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
              kernelreal[i][j] = 0.0;
              kernelimag[i][j] = 0.0;
            }
          }

          for(int theta_idx=0; theta_idx < BSE_ang_absciss_points; theta_idx++){

            theta = BSE_absciss_ang[theta_idx];
            q_p = p*q*(z*std::cos(alpha) + std::sin(psi)*std::cos(theta)*std::sin(alpha));

            k_k = p*p + q*q - 2.0*p*q*(z*std::cos(alpha) + std::sin(psi)*std::cos(theta)*std::sin(alpha));
            k_q = q_q - q_p;
            k_p = q_p - p_p;


            #include "pion_bse_4D.c"

            // for(int dirac_i_idx=0; dirac_i_idx < 2;  dirac_i_idx++){
            //   for(int dirac_j_idx=0; dirac_j_idx < 2; dirac_j_idx++){
            //
            //     if(dirac_i_idx==0 && dirac_j_idx==0){
            //
            //         dressing1sumreal += (kernel[0][0]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).real();
            //         dressing1sumimag += (kernel[0][0]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).imag();
            //
            //           std::complex<double> dressing1 = {dressing1sumreal, dressing1sumimag};
            //
            //           matrix_entry = dressing1;
            //       }
            //
            //       else if(dirac_i_idx==0 && dirac_j_idx==1){
            //
            //           dressing2sumreal += (kernel[0][1]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).real();
            //           dressing2sumimag += (kernel[0][1]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).imag();
            //
            //           std::complex<double> dressing2 = {dressing2sumreal,dressing2sumimag};
            //           matrix_entry = dressing2;
            //           }
            //
            //         else if(dirac_i_idx==1 && dirac_j_idx==0){
            //
            //           dressing3sumreal += (kernel[1][0]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).real();
            //           dressing3sumimag += (kernel[1][0]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).imag();
            //
            //           std::complex<double> dressing3 = {dressing3sumreal,dressing3sumimag};
            //           matrix_entry = dressing3;
            //         }
            //
            //           else{
            //
            //               dressing4sumreal += (kernel[1][1]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).real();
            //               dressing4sumimag += (kernel[1][1]*BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).imag();
            //
            //           std::complex<double> dressing4 = {dressing4sumreal,dressing4sumimag};
            //           matrix_entry = dressing4;
            //         }
            //
            //       alpha_idx = psi_p_idx + p_idx*BSE_ang_absciss_points + dirac_i_idx*BSE_absciss_points*BSE_ang_absciss_points;
            //       beta_idx = psi_idx + q_idx*BSE_ang_absciss_points + dirac_j_idx*BSE_absciss_points*BSE_ang_absciss_points;
            //       mother_temp[alpha_idx][beta_idx] = matrix_entry ;
            //       // #pragma omp critical
            //       // std::cout<<"alpha_idx = "<<alpha_idx << " beta_idx = " << beta_idx<<" dirac_i = "<< dirac_i_idx<< " dirac_j = "<< dirac_j_idx << " matrix_entry = " << matrix_entry<< std::endl;
            //
            //
            //       }
            //     }
              }
            }
          }
        }
      }
      // std::cout<< mother_temp[2040][1200]<< std::endl;
      std::cout<<"MotherMatrix generated" << std::endl;
  return mother_temp;

}


double int_coupled_a(double*** angular_matrix, double* absciss_x, double* weights_w, double* a_vals, double* b_vals, int p_idx, double m_g){
  double s0_a = 0.0;
  double c1,c2; //Prefactor of Integral
  double q, z, yota, k_squared, angularpta;

#pragma omp parallel for private(q) default(none) shared(weights_w, absciss_x, a_vals, b_vals, p_idx, angular_matrix) reduction(+:s0_a)
  for(int q_idx=0;q_idx<absciss_points;q_idx++){

#ifdef MarisTandy

      #ifdef loggrid

      q = exp(0.5*absciss_x[q_idx]);

      s0_a += (weights_w[q_idx] * angular_matrix[0][q_idx][p_idx] * a_vals[q_idx] / (pow(q*a_vals[q_idx],2.0) + pow(b_vals[q_idx],2.0)));

      #else

      c2 = g_squared/(3.0*pow(M_PI,3.0));
      q = absciss_x[q_idx];
      for(int j=0;j<absciss_points;j++){
        z = cos(absciss_ang[j]);
        yota = absciss_points[j];
        k_squared = p*p + q*q - 2.0*p*q*z;
        if(p==0.0){
          return 0.0;
        }
        else{
        s0_a += weights_w[q_idx] *(1.0/(p*p))*( (c2 * q*q*q * a_vals[q_idx]) / ((q*q * pow(a_vals[q_idx],2.0) + pow(b_vals[q_idx],2.0))) ) *
                weights_ang[j] * sin(yota)*sin(yota) * (p*q*z + (2.0/(k_squared)) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z)) *
                (running_coupling_MarisTandy(k_squared,eta) / (k_squared));
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;
        }
      }

      #endif


#else

      c1 = 2.0/(3.0*pow(m_g,2.0)*pow(M_PI,3.0));

      for(int j=0;j<absciss_points;j++){
        z = cos(absciss_ang[j]);
        yota = absciss_ang[j];
        if(p==0.0){
          return 0.0;
        }
        else{
        s0_a += weights_w[q_idx]*((c1 * q*q*q*q*q * a_vals[q_idx])/(p*(q*q*pow(a_vals[q_idx],2)+pow(b_vals[q_idx],2))))*
                (weights_ang[j]*sin(yota)*sin(yota)*z);
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;
        }
      }
#endif
  }
  return s0_a;
}

double int_coupled_b(double*** angular_matrix, double* absciss_x, double* weights_w, double* a_vals, double* b_vals, int p_idx, double m_g){
  double s0_b = 0.0;
  double c1,c2;
  c1 = 2.0/(3.0*pow(m_g,2.0)*pow(M_PI,3.0));
  double q, z, yota, k_squared;

#pragma omp parallel for private(q) default(none) shared(weights_w, absciss_x, a_vals, b_vals, p_idx, angular_matrix) reduction(+:s0_b)
  for(int q_idx=0;q_idx<absciss_points;q_idx++){


#ifdef MarisTandy

      #ifdef loggrid

        q = exp(0.5*absciss_x[q_idx]);

        s0_b += (weights_w[q_idx] * angular_matrix[1][q_idx][p_idx] * b_vals[q_idx] / (pow(q * a_vals[q_idx],2.0) + pow(b_vals[q_idx],2.0)));

      #else

          c2 = (g_squared/pow(M_PI,3.0));

          q = absciss_x[q_idx];

          for(int j=0;j<absciss_points;j++){
            z = cos(absciss_ang[j]);
            yota = absciss_ang[j];
            k_squared = p*p + q*q - 2.0*p*q*z;
            if(p==0.0){
              return 0.0;
            }
            else{
            s0_b += (weights_w[q_idx] * (c2 * q*q*q * b_vals[q_idx] / (q*q * pow(a_vals[q_idx],2) + pow(b_vals[q_idx],2)))) *
                    (weights_ang[j] * (sin(yota)*sin(yota)*(running_coupling_MarisTandy(k_squared, eta)/(k_squared)))) ;
            // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
            // m0B = m_c + c_1*s0_b;
            // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
            // return 0;
            }
          }

      #endif


#else

    for(int j=0;j<absciss_points;j++){

      z = cos(absciss_ang[j]);
      yota = absciss_ang[j];
        // mend=m0B;
        s0_b += weights_w[q_idx]*((c1 * q*q*q*q * b_vals[q_idx])/(q*q*pow(a_vals[q_idx],2)+pow(b_vals[q_idx],2)))*
                (weights_ang[j]*sin(yota)*sin(yota));
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;

      }

#endif

  }
  return s0_b;
}

double** iterate_dressing_functions(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double g_squared, double eta, double mu){

  double*** angular_matrix = initialize_matrix(epsilon,m_c,absciss_x, absciss_ang, weights_ang,g_squared,eta,mu);

  double** init = initialize_dressing_functionAB(1.0,m_c);
  double* a_vals = init[0];
  double* b_vals = init[1];
  double* renorm_constants = init[2];
  double p;

  double new_b[absciss_points];
  double new_a[absciss_points];
  double z2, z2_old, zm, new_siga, new_sigb;
  double renorm_a[absciss_points];
  double renorm_b[absciss_points];

  // double a_end=0.1;
  double a_start=0.0;
  double b_start=0.0;
  double a_end=1.0;
  double b_end=1.0;

  double a_renorm_s=0.0;
  double b_renorm_s=0.0;
  double a_renorm_end=1.0;
  double b_renorm_end=1.0;

  zm = 1.0;
  z2 = 1.0;


  ProgressBar pb(max_iter, "Iterating Dressing Functions A and B");

  for(int j=0;j<=max_iter;j++){
      ++pb;
      if((abs((b_end-b_start)/(b_end+b_start))<epsilon) && (abs((a_end-a_start)/(a_end+a_start))<epsilon)) { //(abs((b_end-b_start)/(b_end+b_start))<epsilon)
        std::cout<<std::endl<<"Success! A and B converged. "<< j << " iterations used. Maximum Iterations were set to "<< max_iter <<std::endl<<std::endl;

        double** dress2d = nullptr;
        dress2d = new double*[3];
        std::swap(dress2d[0],a_vals);
        std::swap(dress2d[1],b_vals);
        std::swap(dress2d[2],renorm_constants);
        return dress2d;
        //
        // for (int i = 0; i<3 ; i++){
        //   free(dress2d[i]);
        // }

      }
      else{

        a_start = new_a[absciss_points-1];
        b_start = new_b[absciss_points-1];
        // ProgressBar ABupdate(absciss_points, "Calculating New A and B");

        for(int p_idx=0;p_idx<absciss_points;p_idx++){ //Integration over q
          // ++ABupdate;

          new_a[p_idx] = renorm_constants[0]*1.0;
          new_b[p_idx] = renorm_constants[0]*renorm_constants[1]*m_c;

          new_a[p_idx] += renorm_constants[0]*renorm_constants[0]*int_coupled_a(angular_matrix, absciss_x, weights_w, a_vals, b_vals, p_idx, m_g);
          new_b[p_idx] += renorm_constants[0]*renorm_constants[0]*int_coupled_b(angular_matrix, absciss_x, weights_w, a_vals, b_vals, p_idx, m_g);

      }

          new_siga = int_coupled_a(angular_matrix, absciss_x, weights_w, a_vals, b_vals, absciss_points, m_g);
          new_sigb = int_coupled_b(angular_matrix, absciss_x, weights_w, a_vals, b_vals, absciss_points, m_g);

          renorm_constants[0] = 1.0/(1.0 + renorm_constants[0]*new_siga);
          renorm_constants[1] = 1.0/renorm_constants[0] - renorm_constants[0]*new_sigb/m_c;



        // ##### UPDATE A and B Values.

        for(int k=0; k < absciss_points; k++){

          std::swap(a_vals[k],new_a[k]);
          std::swap(b_vals[k],new_b[k]);

          // a_vals[k] = new_a[k];
          // b_vals[k] = new_b[k];

        }

        // ##### UPDATE temporarily cross-check values.

      std::swap(a_end,new_a[absciss_points-1]);
      std::swap(b_end,new_b[absciss_points-1]);

    }
  }
  // free(init[0]);
  // free(init[1]);
  // free(init[2]);
  // free(init);
  return 0;
}

std::complex<double>* interpolation_cmplx(std::complex<double> p, double m_c, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x, double* absciss_ang, double* weights_w, double* weights_ang, double eta){
  std::complex<double>* interpolvals = new std::complex<double>[3];
  std::complex<double> sa, sb;

  std::complex<double> k_squared;
  double yota,q,z;
  double constant_a_real, constant_a_imag, constant_b_imag, constant_b_real, sa_imag, sa_real, sb_real, sb_imag;
  sa_real = renorm_constants[0];
  sb_real = renorm_constants[0]*renorm_constants[1]*m_c;
  sa_imag = 0.0;
  sb_imag = 0.0;

  for(int q_idx = 0; q_idx < absciss_points; q_idx++){
      q = std::exp(0.5*absciss_x[q_idx]);

#pragma omp parallel for default(none) shared(p,q,q_idx,absciss_ang,weights_ang,weights_w,renorm_constants,a_vals,b_vals,eta) private(yota, z, k_squared,constant_a_imag,constant_a_real,constant_b_imag,constant_b_real) reduction(+:sa_real, sa_imag,sb_real,sb_imag)

      for(int ang_idx = 0; ang_idx < ang_absciss_points; ang_idx++){
        yota = absciss_ang[ang_idx];
        z = std::cos(yota);
        k_squared = p*p + q*q - 2.0*p*q*z;

        constant_a_real = ((1.0/(p*p)) * (p*q*z + (2.0/(k_squared)) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z)) *
                        (running_coupling_MarisTandy_cmplx(k_squared,eta) / (k_squared))).real();

        constant_a_imag = ((1.0/(p*p)) * (p*q*z + (2.0/(k_squared)) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z)) *
                        (running_coupling_MarisTandy_cmplx(k_squared,eta) / (k_squared))).imag();

        constant_b_real = (running_coupling_MarisTandy_cmplx(k_squared, eta)/(k_squared)).real();

        constant_b_imag = (running_coupling_MarisTandy_cmplx(k_squared, eta)/(k_squared)).imag();


        sa_real += renorm_constants[0]*renorm_constants[0]* weights_w[q_idx]*((4.0/(2.0*3.0*std::pow(M_PI,2.0)) * q*q*q*q * a_vals[q_idx]) /
              ((q*q * std::pow(a_vals[q_idx],2.0) + std::pow(b_vals[q_idx],2.0)))) * constant_a_real *
              weights_ang[ang_idx] * std::sin(yota)*std::sin(yota) ;

        sb_real += renorm_constants[0]*renorm_constants[0]*(weights_w[q_idx] * (4.0/(2.0*std::pow(M_PI,2.0)) * q*q*q*q * b_vals[q_idx] /
              (q*q * std::pow(a_vals[q_idx],2.0) + std::pow(b_vals[q_idx],2.0)))) *
              (weights_ang[ang_idx] * (std::sin(yota)*std::sin(yota)*constant_b_real)) ;

        sa_imag += renorm_constants[0]*renorm_constants[0]* weights_w[q_idx]*((4.0/(2.0*3.0*std::pow(M_PI,2.0)) * q*q*q*q * a_vals[q_idx]) /
              ((q*q * std::pow(a_vals[q_idx],2.0) + std::pow(b_vals[q_idx],2.0)))) * constant_a_imag *
              weights_ang[ang_idx] * std::sin(yota)*std::sin(yota) ;

        sb_imag += renorm_constants[0]*renorm_constants[0]*(weights_w[q_idx] * (4.0/(2.0*std::pow(M_PI,2.0)) * q*q*q*q * b_vals[q_idx] /
              (q*q * std::pow(a_vals[q_idx],2.0) + std::pow(b_vals[q_idx],2.0)))) *
              (weights_ang[ang_idx] * (std::sin(yota)*std::sin(yota)*constant_b_imag)) ;

        // fileoutD<< q << "\t" << z << "\t" << sa_real << "\t" << sa_imag << "\t" << sb_real <<"\t" << sb_imag << std::endl;

      }
  }

  sa = {sa_real,sa_imag};
  sb = {sb_real,sb_imag};

  interpolvals[0] = sa;
  interpolvals[1] = sb;
  interpolvals[2] = interpolvals[1]/interpolvals[0];

  return interpolvals;
}

// ##### FREE STUFF ##### //

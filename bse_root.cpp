#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "definitions.h"
#include "help_numerics.h"
#include "Dyson_test.h"
#include "progressbar.hpp"
#include <complex>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>  // Also includes <MatOp/DenseGenMatProd.h>

double regulaFalsi(double low_mass, double high_mass, double epsilon, std::complex<double>* a_corners, std::complex<double>* b_corners, double* BSE_weights_ang, double* BSE_absciss_x,
  double* BSE_absciss_ang, std::complex<double>***** theta_matrix, double* q_vec, double* z_vec, double* x_corners, std::complex<double>*** y_corner)
{
      double f_low = bse_root_eigenlib(low_mass,a_corners,  b_corners, BSE_weights_ang, BSE_absciss_x, BSE_absciss_ang,theta_matrix, q_vec, z_vec, x_corners, y_corner);
      double f_high = bse_root_eigenlib(high_mass,a_corners,  b_corners, BSE_weights_ang, BSE_absciss_x, BSE_absciss_ang, theta_matrix, q_vec, z_vec, x_corners, y_corner);

      double x_low, x_high, delta_x, delta_mass, c, possibleroot;

    if ( f_low * f_high >= 0.0)
    {
        std::cout << "Variables low_mass and high_mass not assumed correctly. Try again\n";
        return 0;
    }

    if(f_low < 0.0){
      x_low = low_mass;
      x_high = high_mass;
    }
    else{
      x_low = high_mass;
      x_high = low_mass;
    }

    delta_x = high_mass - low_mass;

    for (int i=0; i < max_iter; i++)
    {
        // Find the point that touches x axis
        c = ((x_low*f_high - x_high*f_low )/ (f_high - f_low));
        possibleroot = bse_root_eigenlib(c, a_corners,  b_corners, BSE_weights_ang, BSE_absciss_x, BSE_absciss_ang, theta_matrix, q_vec, z_vec, x_corners, y_corner);

        // Check if the above found point is root
        if(possibleroot*f_high < 0.0){
          delta_mass = x_low - c;
          x_low = c;
          f_low = possibleroot;
        }
        else{
          delta_mass = x_high - c;
          x_high = c;
          f_high = possibleroot;
        }
        delta_x = x_high - x_low;

        // std::cout<<std::endl<<possibleroot<<std::endl;
        // Decide the side to repeat the steps
        if(abs(possibleroot) < epsilon){
          std::cout<< "Difference to iterated 'exact' value: "<< abs(possibleroot) - epsilon << std::endl;

          // std::cout << "The value of root is : " << c << std::endl;
          return c;
          break;
        }
    }
    std::cout<< "Reached maximum Iterations for regulaFalsitest! " <<std::endl;
    return 0;
}


double regulaFalsitest(double low_mass, double high_mass, double epsilon)
{
      double f_low = low_mass*low_mass-4.0*low_mass-10.0;
      double f_high = high_mass*high_mass-4.0*high_mass-10.0;
      double x_low, x_high, delta_x, delta_mass, c, possibleroot;

    if ( f_low * f_high >= 0.0)
    {
        std::cout << "Variables low_mass and high_mass not assumed correctly. Try again\n";
        return 0;
    }

    if(f_low < 0.0){
      x_low = low_mass;
      x_high = high_mass;
    }
    else{
      x_low = high_mass;
      x_high = low_mass;
    }

    delta_x = high_mass - low_mass;

    for (int i=0; i < max_iter; i++)
    {
        // Find the point that touches x axis
        c = ((x_low*f_high - x_high*f_low )/ (f_high - f_low));
        possibleroot = c*c-4.0*c-10.0;

        // Check if the above found point is root
        if(possibleroot*f_high < 0.0){
          delta_mass = x_low - c;
          x_low = c;
          f_low = possibleroot;
        }
        else{
          delta_mass = x_high - c;
          x_high = c;
          f_high = possibleroot;
        }
        delta_x = x_high - x_low;

        std::cout<<std::endl<<possibleroot<<std::endl;
        // Decide the side to repeat the steps
        if(abs(possibleroot) < epsilon){

          std::cout << "The value of root is : " << c << std::endl;
          return c;
          break;
        }
    }
    std::cout<< "Reached maximum Iterations for regulaFalsitest! " <<std::endl;
    return 0;
}


struct quadratic_params
  {
    double m_c,eta;
    double* renorm_constants,a_vals,b_vals,absciss_x,absciss_ang,weights_w,weights_ang;
  };

// double bsefunction(double mass, void* params){
//   struct quadratic_params *p = (struct quadratic_params *) params;
//
//   double m_c = p-> m_c;
//   double eta = p-> eta;
//   double* renorm_constants = p-> renorm_constants;
//   double* a_vals = p-> a_vals;
//   double* b_vals = p-> b_vals;
//   double* absciss_x = p-> absciss_x;
//   double* absciss_ang = p-> absciss_ang;
//   double* weights_w = p-> weights_w;
//   double* weights_ang = p-> weights_ang;
//
//   return bse_root(mass, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x,
//      absciss_ang,  weights_w,  weights_ang,  eta);
// }
//
// double findRoot(double low_mass, double high_mass, double m_c, double eta, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x,
//   double* absciss_ang, double* weights_w, double* weights_ang){
//
// int status;
//   int iter = 0;
//   const gsl_root_fsolver_type *T;
//   gsl_root_fsolver *s;
//   double r = 0, r_expected = sqrt (5.0);
//   double x_lo = low_mass, x_hi = high_mass;
//   gsl_function F;
//   struct quadratic_params params = {m_c,eta,renorm_constants,a_vals,b_vals,absciss_x,absciss_ang,weights_w, weights_ang};
//
//   F.function = &bsefunction;
//   F.params = &params;
//
//   T = gsl_root_fsolver_brent;
//   s = gsl_root_fsolver_alloc (T);
//   gsl_root_fsolver_set (s, &F, x_lo, x_hi);
//
//   printf ("using %s method\n",
//           gsl_root_fsolver_name (s));
//
//   printf ("%5s [%9s, %9s] %9s %10s %9s\n",
//           "iter", "lower", "upper", "root",
//           "err", "err(est)");
//
//   do
//     {
//       iter++;
//       status = gsl_root_fsolver_iterate (s);
//       r = gsl_root_fsolver_root (s);
//       x_lo = gsl_root_fsolver_x_lower (s);
//       x_hi = gsl_root_fsolver_x_upper (s);
//       status = gsl_root_test_interval (x_lo, x_hi,
//                                        0, 0.001);
//
//       if (status == GSL_SUCCESS)
//         printf ("Converged:\n");
//
//       printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//               iter, x_lo, x_hi,
//               r, r - r_expected,
//               x_hi - x_lo);
//     }
//   while (status == GSL_CONTINUE && iter < max_iter);
//
//   gsl_root_fsolver_free (s);
//
//   return status;
// }

// using namespace Spectra;
double bse_root_eigenlib(double pionmass, std::complex<double>* a_corners, std::complex<double>* b_corners, double* BSE_weights_ang, double* BSE_absciss_x,
  double* BSE_absciss_ang, std::complex<double>***** theta_matrix, double* q_vec, double* z_vec, double* x_corners, std::complex<double>*** y_corner){


    // std::cout<<std::endl<< "Generating Mother Matrix" << std::endl;
    std::complex<double>** mother = initialize_mother_matrix(pionmass, a_corners, b_corners, BSE_weights_ang, BSE_absciss_x, BSE_absciss_ang, theta_matrix, q_vec, z_vec, x_corners, y_corner);
    int grid = BSE_absciss_points*BSE_ang_absciss_points*4;
    double root;

    // for(int i=0; i<grid; i++){
    //   for(int j=0; j<grid; j++){
    //     std::cout<<"alpha_idx = "<< i << " beta_idx = " << j << mother[i][j] << std::endl;
    //   }
    // }
    Eigen::MatrixXd Ematrix(grid,grid);
    // Ematrix = mother.real();


    for(int i=0; i<grid; i++){
      for(int j=0; j<grid; j++){
        Ematrix(i,j) = (mother[i][j]).real();
      }
    }
    Spectra::DenseGenMatProd<double> op(Ematrix);
    // std::cout<<" HERE OK " << std::endl;

    Spectra::GenEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::DenseGenMatProd<double> > eigs(&op,1,6);

    eigs.init();
    int nconv = eigs.compute();

    Eigen::VectorXcd evalues;
    if(eigs.info() == Spectra::SUCCESSFUL)
        evalues = eigs.eigenvalues();
        root = evalues[0].real();
    std::cout << "Eigenvalues found:\n" << root << std::endl;

    return root - 1.0;

    // Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    // ces.compute(Ematrix);
    // for(int j=0;j<absciss_points;j++){
    //   // ++sd;
    //    std::cout << "The eigenvalues of A are:" << std::endl << ces.eigenvalues()[j] << std::endl;
    // }
    // std::cout << "The matrix of eigenvectors, V, is:" << std::endl << ces.eigenvectors().col(absciss_points -1).real() << std::endl << std::endl;
    //
    // double root = ces.eigenvalues()[absciss_points - 1].real();
    // std::cout<<root<<std::endl;
    //
    // std::ofstream  pionmassfile;
    // pionmassfile.open("Data/PionmassandE(pP).dat");
    // pionmassfile << "# Parameters used: " << "mc(GeV): "<< m_c<<" LAMBDA(GeV) in UV-Cuttoff in log(LAMBDA*LAMBDA): "<< LAMBDA << "LAMBDA_MIN(GeV) in IR-Cuttoff in log(LAMBDA_MIN*LAMBDA_MIN): "<< LAMBDA_MIN<< " gamma_m: "<< gamma_fun(N_C, N_F) <<std::endl;
    // pionmassfile << "# mu(GeV): "<< mu_renorm << " Lambda_QCD(GeV): "<< Lambda_QCD << " Lambda_t(GeV): "<< Lambda_t << " Lambda_0 " <<Lambda_0<<std::endl;
    // pionmassfile << "# q-abscissae used: "<< absciss_points <<" ang_abscissae used: "<< ang_absciss_points <<std::endl;
    // pionmassfile << "# z2 is: " << renorm_constants[0] << " zm is: " << renorm_constants[1]<<std::endl;
    // pionmassfile << "# Pionmass is : " << pionmass<<std::endl;
    // pionmassfile << "# p^2"<< " "<< "E(p,P)"<<std::endl;
    // for(int j=0;j<absciss_points;j++){
    //   // ++sd;
    //   pionmassfile<< exp(absciss_x[j]) << " " << ces.eigenvectors().col(absciss_points-1).real()[j]<< std::endl;
    // }
    // pionmassfile.close ();

    // return 0;

  }

double bse_root(double pionmass, std::complex<double>* a_corners, std::complex<double>* b_corners, double* BSE_weights_ang, double* BSE_absciss_x,
  double* BSE_absciss_ang, std::complex<double>***** theta_matrix, double* q_vec, double* z_vec, double* x_corners, std::complex<double>*** y_corner){

std::cout<<std::endl<< "Generating Mother Matrix" << std::endl;
std::complex<double>** mother = initialize_mother_matrix(pionmass, a_corners, b_corners, BSE_weights_ang, BSE_absciss_x, BSE_absciss_ang, theta_matrix, q_vec, z_vec, x_corners, y_corner);
int grid = BSE_absciss_points*BSE_ang_absciss_points*4;
// double* gsl_mother_temp = nullptr;
// gsl_mother_temp = new double[absciss_points*absciss_points];

// for(int i = 0; i< absciss_points; i ++){
//   for(int j = 0; j< absciss_points; j++){
//     std::cout<< mother[i][j].real() << " " << mother[i][j].imag() << std::endl;
//   }
// }



// ##### WRITE FILE ######
// std::ofstream fileout3;
// fileout3.open("Data/MotherMatrix_temp.dat");
// for(int i = 0; i< absciss_points; i ++){
//   for(int j = 0; j< absciss_points; j++){
//     fileout3 << mother[i][j].real() << std::endl;
//   }
// }
// fileout3.close();
//
// // ##### READ FILE #####
// std::ifstream input("Data/MotherMatrix_temp.dat");
//  if (!input) {
//    std::cout << "Cannot open file.\n";
//    return 0;
//  }
//  for (int i = 0; i < absciss_points*absciss_points; i++) {
//    input >> gsl_mother_temp[i];
//  }
//  input.close();
// #####  #####

// ##### SOLVE EIGENVALUE PROBLEM WITH GSL LIBRARY #####
// gsl_matrix_view gsl_mother
//   = gsl_matrix_view_array (gsl_mother_temp, absciss_points, absciss_points);

gsl_vector_complex *alpha = gsl_vector_complex_alloc (grid);
gsl_vector* beta = gsl_vector_alloc(grid);
gsl_matrix_complex *evec = gsl_matrix_complex_alloc (grid, grid);

gsl_eigen_genv_workspace * mother_workspace =
  gsl_eigen_genv_alloc (grid);

gsl_matrix* gsl_mother_temp = gsl_matrix_alloc(grid,grid);
gsl_matrix* one = gsl_matrix_alloc(grid,grid);

#pragma omp parallel for default(none) shared(gsl_mother_temp,grid,one,mother)
for(int i=0; i<grid; i++){
  for(int j=0; j<grid; j++){
    gsl_matrix_set(gsl_mother_temp,i,j,mother[i][j].real());
    gsl_matrix_set(one,i,j,0.0);
    if(i==j){
      gsl_matrix_set(one,i,j,1.0);
    }
  }
}

// gsl_eigen_nonsymmv_params(0, 1, mother_workspace);
gsl_eigen_genv (gsl_mother_temp,one, alpha,beta,evec, mother_workspace);

gsl_eigen_genv_free (mother_workspace);

gsl_eigen_genv_sort (alpha,beta,evec,
                         GSL_EIGEN_SORT_ABS_DESC);

double alpha0 = GSL_REAL(gsl_vector_complex_get(alpha,0));
double beta0 = gsl_vector_get(beta,0);
double root = alpha0/beta0;


std::cout<< "first eigenvalue of mother-matrix is: "<<root<<std::endl;
// std::cout<< "the others are: " <<std::endl;
// for (int j = 0; j < absciss_points; ++j)
//           {
//             std::cout<< j << " " << GSL_REAL(gsl_vector_complex_get(eval,j)) << " + i* "<< GSL_IMAG(gsl_vector_complex_get(eval,j))<<std::endl;
//           }
//           std::cout<<std::endl;



gsl_vector_complex_view evec_0
           = gsl_matrix_complex_column (evec, 0);

// std::cout<<"geile eigenvectors:"<<std::endl;
//    for(int j = 0; j<absciss_points; j++){
//      std::cout<< GSL_REAL(gsl_matrix_complex_get(evec,j,0))<<std::endl;
//    }

std::cout<<"vector again"<<std::endl;
for (int j = 0; j < grid; ++j)
          {
            gsl_complex z =
              gsl_vector_complex_get(&evec_0.vector, j);
            std::cout<<GSL_REAL(z)<<" + i* "<< GSL_IMAG(z)<<std::endl;
          }

return root - 1.0;



gsl_vector_complex_free(alpha);
gsl_vector_free(beta);
gsl_matrix_complex_free(evec);
gsl_matrix_free(gsl_mother_temp);

// return 0;
}

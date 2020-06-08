
    w[1]=1.0/( - P_P*p_p + pow(P_p,2));
    w[2]=pow(k_k,-1);
    w[3]=k_P;
    w[4]=k_q;
    w[5]=p_p;
    w[6]=k_p;
    w[7]=P_p;
    w[8]=P_P;
    w[9]=P_q;
    w[10]=q_p;
    w[11]=1.0/( - P_P*pow(P_p,2)*p_p + pow(P_p,4));
    w[12]=1.0/( - 1.6E+1*P_P*p_p + 1.6E+1*pow(P_p,2));
    w[13]=q_q;
   w[14]=sigmasm*sigmavp;
   w[15]=sigmasp*sigmavm;
   w[16]=w[14] - w[15];
   w[17]=w[16]*w[9];
   w[18]=w[16]*eta;
   w[18]=w[18] + w[15];
   w[19]=w[18]*w[8];
   w[17]=w[17] + w[19];
   w[19]= - w[5]*w[17];
   w[20]=w[10]*w[16];
   w[21]=w[7]*w[18];
   w[20]=w[20] + w[21];
   w[20]=w[7]*w[20];
   w[21]=w[16]*w[4];
   w[22]=w[18]*w[3];
   w[21]=w[21] + w[22];
   w[22]=w[5]*w[3];
   w[23]=w[7]*w[6];
   w[22]=w[22] - w[23];
   w[23]=2.E+0*w[2];
   w[22]=w[22]*w[23];
   w[24]= - w[21]*w[22];
   w[19]=w[24] + w[19] + w[20];
   w[19]=w[1]*w[19];
   w[20]=w[10]*w[8];
   w[24]=w[7]*w[9];
   w[25]= - w[20] + w[24];
   w[25]=w[7]*w[16]*w[25];
   w[26]=w[7]*w[3];
   w[27]= - w[6]*w[8];
   w[27]=w[27] + w[26];
   w[28]=w[23]*w[7];
   w[21]=w[28]*w[21]*w[27];
   w[21]=w[25] + w[21];
   w[21]=w[11]*w[21];
   w[25]=w[26]*w[4];
   w[27]=w[10]*pow(w[3],2);
   w[25]=w[25] - w[27];
   w[27]=sigmavp*sigmavm;
   w[29]= - w[27]*w[25];
   w[30]=w[4]*w[8];
   w[31]=w[27]*w[30];
   w[32]=w[27]*w[9];
   w[33]= - w[3]*w[32];
   w[31]=w[31] + w[33];
   w[31]=w[6]*w[31];
   w[29]=w[31] + w[29];
   w[29]=w[29]*w[23];
   w[31]= - w[27]*w[20];
   w[33]=w[7]*w[32];
   w[29]=w[29] + w[31] + w[33];
   w[31]=8.E+0*w[12];
   w[29]=w[29]*w[31];
   w[17]=3.E+0*w[17];
   w[33]=eta - 1.E+0;
   w[33]=w[27]*w[33];
   w[34]=w[33]*w[8]*eta;
   w[35]=sigmasm*sigmasp;
   w[36]=w[27]*w[13];
   w[37]=w[35] - w[36];
   w[38]=w[34] + w[37];
   w[39]=w[8]*w[38];
   w[40]=2.E+0*eta;
   w[41]=w[40] - 1.E+0;
   w[41]=w[27]*w[41];
   w[42]=w[41]*w[8];
   w[32]=w[42] + 2.E+0*w[32];
   w[43]=w[9]*w[32];
   w[39]=w[39] + w[43];
   w[39]=w[5]*w[39];
   w[43]= - w[7]*w[38];
   w[44]= - w[10]*w[32];
   w[43]=w[44] + w[43];
   w[43]=w[7]*w[43];
   w[38]=w[38]*w[3];
   w[44]=w[32]*w[4];
   w[44]=w[44] + w[38];
   w[45]=w[44]*w[22];
   w[39]=w[45] + w[39] + w[43];
   w[39]=w[1]*w[39];
   w[43]=w[9]*w[8];
   w[45]=w[43]*w[27];
   w[46]=pow(w[8],2);
   w[47]=w[46]*w[41];
   w[45]=w[47] + 2.E+0*w[45];
   w[47]=w[4]*w[45];
   w[38]=w[8]*w[38];
   w[38]=w[47] + w[38];
   w[38]=w[6]*w[38];
   w[44]= - w[44]*w[26];
   w[38]=w[38] + w[44];
   w[38]=w[38]*w[28];
   w[44]=w[10]*w[45];
   w[32]= - w[32]*w[24];
   w[32]=w[44] + w[32];
   w[32]=w[7]*w[32];
   w[32]=w[32] + w[38];
   w[32]=w[11]*w[32];
   w[14]=w[14] + w[15];
   w[38]=w[14]*w[25];
   w[30]= - w[14]*w[30];
   w[44]=w[14]*w[9];
   w[45]=w[3]*w[44];
   w[30]=w[30] + w[45];
   w[30]=w[6]*w[30];
   w[30]=w[30] + w[38];
   w[30]=w[30]*w[23];
   w[20]=w[14]*w[20];
   w[38]= - w[7]*w[44];
   w[20]=w[30] + w[20] + w[38];
   w[20]=w[20]*w[31];
   w[16]=w[13]*w[16];
   w[18]=w[9]*w[18];
   w[16]=w[16] + w[18];
   w[16]=3.E+0*w[9]*w[16];
   w[18]=w[35] + w[36];
   w[30]=w[34] - w[18];
   w[35]=w[10]*w[30]*w[9];
   w[33]=w[33]*w[40];
   w[38]=w[33]*w[9];
   w[40]=w[41]*w[13];
   w[38]=w[38] + w[40];
   w[40]= - w[38]*w[24];
   w[35]=w[35] + w[40];
   w[35]=w[7]*w[35];
   w[40]=w[3]*w[9];
   w[38]=w[38]*w[40];
   w[45]=w[30]*w[4];
   w[47]=w[45]*w[9];
   w[38]=w[38] - w[47];
   w[47]=w[38]*w[22];
   w[18]=w[34] + w[18];
   w[48]=w[9]*w[18];
   w[42]=w[42]*w[13];
   w[48]=w[42] + w[48];
   w[48]=w[5]*w[9]*w[48];
   w[35]=w[47] + w[48] + w[35];
   w[35]=w[1]*w[35];
   w[45]= - w[43]*w[45];
   w[33]=w[43]*w[33];
   w[33]=w[42] + w[33];
   w[33]=w[33]*w[40];
   w[33]=w[45] + w[33];
   w[33]=w[6]*w[33];
   w[38]= - w[38]*w[26];
   w[33]=w[33] + w[38];
   w[33]=w[33]*w[28];
   w[38]=pow(w[9],2);
   w[42]=w[38]*w[7];
   w[45]= - w[10]*w[43];
   w[42]=w[45] + w[42];
   w[30]=w[7]*w[30]*w[42];
   w[30]=w[30] + w[33];
   w[30]=w[11]*w[30];
   w[33]=w[14]*eta;
   w[15]=w[33] - w[15];
   w[33]=w[15]*w[9];
   w[42]= - w[33]*w[25];
   w[45]=w[15]*w[43];
   w[47]=w[4]*w[45];
   w[48]=w[38]*w[15];
   w[49]= - w[3]*w[48];
   w[47]=w[47] + w[49];
   w[47]=w[6]*w[47];
   w[42]=w[47] + w[42];
   w[42]=w[42]*w[23];
   w[47]= - w[10]*w[45];
   w[48]=w[7]*w[48];
   w[42]=w[42] + w[47] + w[48];
   w[31]=w[42]*w[31];
   w[27]=w[38]*w[27];
   w[36]= - w[8]*w[36];
   w[27]=w[36] + w[27];
   w[27]=6.E+0*w[27];
   w[36]=w[15]*w[8];
   w[36]=w[36] + w[44];
   w[42]= - w[10]*w[36];
   w[44]=w[14]*w[13];
   w[33]=w[33] + w[44];
   w[47]=w[7]*w[33];
   w[42]=w[42] + w[47];
   w[42]=w[7]*w[42];
   w[33]=w[33]*w[3];
   w[47]=w[36]*w[4];
   w[33]=w[33] - w[47];
   w[22]= - w[33]*w[22];
   w[44]=w[44]*w[8];
   w[38]=w[14]*w[38];
   w[38]= - w[44] + w[38];
   w[38]=w[5]*w[38];
   w[22]=w[22] + w[38] + w[42];
   w[22]=2.E+0*w[1]*w[22];
   w[14]=w[43]*w[14];
   w[15]=w[46]*w[15];
   w[14]=w[14] + w[15];
   w[15]=w[4]*w[14];
   w[38]= - w[44] - w[45];
   w[38]=w[3]*w[38];
   w[15]=w[15] + w[38];
   w[15]=w[6]*w[15];
   w[26]=w[33]*w[26];
   w[15]=w[15] + w[26];
   w[15]=w[15]*w[28];
   w[14]=w[10]*w[14];
   w[26]= - w[36]*w[24];
   w[14]=w[14] + w[26];
   w[14]=w[7]*w[14];
   w[14]=w[14] + w[15];
   w[14]=2.E+0*w[11]*w[14];
   w[15]=w[34] - w[37];
   w[26]=w[15]*w[8];
   w[28]=w[41]*w[43];
   w[26]=w[28] + w[26];
   w[28]=w[4]*w[26];
   w[33]=w[41]*w[9];
   w[15]=w[15] + w[33];
   w[34]= - w[15]*w[40];
   w[28]=w[28] + w[34];
   w[28]=w[6]*w[28];
   w[25]= - w[15]*w[25];
   w[25]=w[28] + w[25];
   w[23]=w[25]*w[23];
   w[25]= - w[10]*w[26];
   w[15]=w[15]*w[24];
   w[15]=w[23] + w[25] + w[15];
   w[15]=1.6E+1*w[12]*w[15];
   w[18]=w[33] + w[18];
   w[18]=3.E+0*w[18];

for(int dirac_i_idx=0; dirac_i_idx < 4; dirac_i_idx++){
 for(int dirac_j_idx=0; dirac_j_idx < 4; dirac_j_idx++){


	 if(dirac_i_idx==0 && dirac_j_idx==0){

	 kernelreal[0][0] += (w[18] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[0][0] += (w[18] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing00 ={kernelreal[0][0],kernelimag[0][0] };
	 matrix_entry = dressing00;
		 }


	 if(dirac_i_idx==0 && dirac_j_idx==1){

	 kernelreal[0][1] += (w[17] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[0][1] += (w[17] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing01 ={kernelreal[0][1],kernelimag[0][1] };
	 matrix_entry = dressing01;
		 }


	 if(dirac_i_idx==0 && dirac_j_idx==2){

	 kernelreal[0][2] += (w[16] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[0][2] += (w[16] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing02 ={kernelreal[0][2],kernelimag[0][2] };
	 matrix_entry = dressing02;
		 }


	 if(dirac_i_idx==0 && dirac_j_idx==3){

	 kernelreal[0][3] += (w[27] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[0][3] += (w[27] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing03 ={kernelreal[0][3],kernelimag[0][3] };
	 matrix_entry = dressing03;
		 }


	 if(dirac_i_idx==1 && dirac_j_idx==0){

	 kernelreal[1][0] += (w[19] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[1][0] += (w[19] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing10 ={kernelreal[1][0],kernelimag[1][0] };
	 matrix_entry = dressing10;
		 }


	 if(dirac_i_idx==1 && dirac_j_idx==1){

	 kernelreal[1][1] += (w[39] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[1][1] += (w[39] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing11 ={kernelreal[1][1],kernelimag[1][1] };
	 matrix_entry = dressing11;
		 }


	 if(dirac_i_idx==1 && dirac_j_idx==2){

	 kernelreal[1][2] += (w[35] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[1][2] += (w[35] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing12 ={kernelreal[1][2],kernelimag[1][2] };
	 matrix_entry = dressing12;
		 }


	 if(dirac_i_idx==1 && dirac_j_idx==3){

	 kernelreal[1][3] += (w[22] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[1][3] += (w[22] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing13 ={kernelreal[1][3],kernelimag[1][3] };
	 matrix_entry = dressing13;
		 }


	 if(dirac_i_idx==2 && dirac_j_idx==0){

	 kernelreal[2][0] += (w[21] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[2][0] += (w[21] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing20 ={kernelreal[2][0],kernelimag[2][0] };
	 matrix_entry = dressing20;
		 }


	 if(dirac_i_idx==2 && dirac_j_idx==1){

	 kernelreal[2][1] += (w[32] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[2][1] += (w[32] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing21 ={kernelreal[2][1],kernelimag[2][1] };
	 matrix_entry = dressing21;
		 }


	 if(dirac_i_idx==2 && dirac_j_idx==2){

	 kernelreal[2][2] += (w[30] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[2][2] += (w[30] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing22 ={kernelreal[2][2],kernelimag[2][2] };
	 matrix_entry = dressing22;
		 }


	 if(dirac_i_idx==2 && dirac_j_idx==3){

	 kernelreal[2][3] += (w[14] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[2][3] += (w[14] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing23 ={kernelreal[2][3],kernelimag[2][3] };
	 matrix_entry = dressing23;
		 }


	 if(dirac_i_idx==3 && dirac_j_idx==0){

	 kernelreal[3][0] += (w[29] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[3][0] += (w[29] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing30 ={kernelreal[3][0],kernelimag[3][0] };
	 matrix_entry = dressing30;
		 }


	 if(dirac_i_idx==3 && dirac_j_idx==1){

	 kernelreal[3][1] += (w[20] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[3][1] += (w[20] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing31 ={kernelreal[3][1],kernelimag[3][1] };
	 matrix_entry = dressing31;
		 }


	 if(dirac_i_idx==3 && dirac_j_idx==2){

	 kernelreal[3][2] += (w[31] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[3][2] += (w[31] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing32 ={kernelreal[3][2],kernelimag[3][2] };
	 matrix_entry = dressing32;
		 }


	 if(dirac_i_idx==3 && dirac_j_idx==3){

	 kernelreal[3][3] += (w[15] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).real();

	 kernelimag[3][3] += (w[15] *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][
psi_idx][p_idx][q_idx]).imag();

	 std::complex<double> dressing33 ={kernelreal[3][3],kernelimag[3][3] };
	 matrix_entry = dressing33;
		 }


	 alpha_idx = psi_p_idx + p_idx*BSE_ang_absciss_points + dirac_i_idx
 *BSE_absciss_points*BSE_ang_absciss_points;
	 beta_idx = psi_idx + q_idx*BSE_ang_absciss_points + dirac_j_idx
 *BSE_absciss_points*BSE_ang_absciss_points;
	 mother_temp[alpha_idx][beta_idx] = matrix_entry;
	 }
}

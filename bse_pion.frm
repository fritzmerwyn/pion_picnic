Off statistics;

ExtraSymbols array,w;

Symbols sigmavp,sigmasp,sigmavm,sigmasm,eta,x;
Vectors k,P,q,p;
Indices mu,nu,i,j;
Functions Kernel,K,pBasis,pProj,Sp,Sm;

Local Kall = sum_(i,0,3,sum_(j,0,3,Kernel(i,j) * x^(i+4*j)));

Local Kernelino = (d_(mu,nu)-(k(mu)*k(nu)/(k.k)))*
                pProj(0) * g_(0,mu) * Sp * pBasis(0) * Sm * g_(0,nu);

id Kernel(i?,j?) = -1*(d_(mu,nu)-(k(mu)*k(nu)/(k.k)))*pProj(i) * g_(0,mu) * Sp * pBasis(j) * Sm * g_(0,nu);

id Sp = sigmasp-i_*sigmavp*(eta*g_(0,P)+g_(0,q));
id Sm = sigmasm-i_*sigmavm*((-1+eta)*g_(0,P)+g_(0,q));

.sort

id pBasis(0) = g5_(0);
id pBasis(1) = (-i_)*g5_(0)*g_(0,P);
id pBasis(2) = (-i_)*P.q*g5_(0)*g_(0,q);
id pBasis(3) = g5_(0)*(g_(0,q)*g_(0,P) - g_(0,P)*g_(0,q));

id pProj(0) = g5_(0)/4;
id pProj(1) = ((-i_/4)*P.p*g5_(0)*g_(0,p))/((P.p)^2-p.p*P.P)+((i_/4)*p.p*g5_(0)*g_(0,P))/((P.p)^2-p.p*P.P);
id pProj(2) = ((i_/4)*P.P*g5_(0)*g_(0,p))/(P.p*((P.p)^2-p.p*P.P))-((i_/4)*g5_(0)*g_(0,P))/((P.p)^2-p.p*P.P);
id pProj(3) = g5_(0)*(g_(0,p)*g_(0,P) - g_(0,P)*g_(0,p))/(16*((P.p)^2-p.p*P.P));

*id k = q-p;

contract;

format float;
format C;
format O4;

trace4,0;
*Print Kernelino;
.sort

b x;                                  * otherwise the optimization will include x's into the array variables *

.sort                                 * prepare Kall for the optimization *

#optimize Kall

b x;                                  * re-bracketing after optimization, because optimization changes the Kall-expression. bracketing allows us to separate the powers of x in the next step *

.sort                                 * necessary, because of the order of commands (bracketing after locals) *

#do i=0,3
#do j=0,3
l Kall'i''j' = Kall[x^{'i'+ 4*'j'}];   * collect coefficients of powers of x *
#enddo
#enddo

.sort

#write <pion_bse_2D.c> "%O"


#write <pion_bse_2D.c> "for(int dirac_i_idx=0; dirac_i_idx < 4; dirac_i_idx++){"
#write <pion_bse_2D.c> " for(int dirac_j_idx=0; dirac_j_idx < 4; dirac_j_idx++){"
#write <pion_bse_2D.c> "\n"

#do i=0,3
#do j=0,3

#write <pion_bse_2D.c> "\t if(dirac_i_idx=='i' && dirac_j_idx=='j'){\n"
#write <pion_bse_2D.c> "\t kernelreal['i']['j'] += (%E *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).real(); \n" , Kall'i''j'
#write <pion_bse_2D.c> "\t kernelimag['i']['j'] += (%E *BSE_weights_ang[theta_idx]*std::sin(theta)*theta_matrix[1][theta_idx][psi_idx][p_idx][q_idx]).imag(); \n" , Kall'i''j'
#write <pion_bse_2D.c> "\t std::complex<double> dressing%s ={kernelreal['i']['j'],kernelimag['i']['j'] };",'i''j'
#write <pion_bse_2D.c> "\t matrix_entry = dressing%s;", 'i''j'
#write <pion_bse_2D.c> "\t\t } \n\n"

#enddo
#enddo

#write <pion_bse_2D.c> "\t alpha_idx = psi_p_idx + p_idx*BSE_ang_absciss_points + dirac_i_idx \n *BSE_absciss_points*BSE_ang_absciss_points;"
#write <pion_bse_2D.c> "\t beta_idx = psi_idx + q_idx*BSE_ang_absciss_points + dirac_j_idx \n *BSE_absciss_points*BSE_ang_absciss_points;"
#write <pion_bse_2D.c> "\t mother_temp[alpha_idx][beta_idx] = matrix_entry;"

#write <pion_bse_2D.c> "\t }"
#write <pion_bse_2D.c> "}"


.end

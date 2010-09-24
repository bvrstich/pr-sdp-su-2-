#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int TPM::counter = 0;

int ***TPM::t2s;
int ***TPM::s2t;

double **TPM::_6j;

/**
 * standard constructor for a spinsymmetrical tp matrix: constructs BlockMatrix object with 2 blocks, for S = 0 or 1,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and tp basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
TPM::TPM(int M,int N) : BlockMatrix(2) {

   this->N = N;
   this->M = M;

   //set the dimension and degeneracy of the two blocks:
   this->setMatrixDim(0,M*(M + 2)/8,1);
   this->setMatrixDim(1,M*(M - 2)/8,3);

   if(counter == 0){

      //allocatie van s2t
      s2t = new int ** [2];

      for(int i = 0;i < 2;++i){

         s2t[i] = new int * [M/2];
         s2t[i][0] = new int [M*M/4];

         for(int j = 1;j < M/2;++j)
            s2t[i][j] = s2t[i][j - 1] + M/2;

      }

      //allocatie van t2s
      t2s = new int ** [2];

      for(int i = 0;i < 2;++i){

         t2s[i] = new int * [this->gdim(i)];

         for(int j = 0;j < this->gdim(i);++j)
            t2s[i][j] = new int [2];

      }

      //initialisatie van de arrays
      int teller = 0;

      //symmetrical array: S = 0
      for(int a = 0;a < M/2;++a)
         for(int b = a;b < M/2;++b){

            s2t[0][a][b] = teller;

            t2s[0][teller][0] = a;
            t2s[0][teller][1] = b;

            ++teller;

         }

      //watch it!
      teller = 0;

      //antisymmetrical array: S = 1
      for(int a = 0;a < M/2;++a)
         for(int b = a + 1;b < M/2;++b){

            s2t[1][a][b] = teller;

            t2s[1][teller][0] = a;
            t2s[1][teller][1] = b;

            ++teller;

         }

      //symmetrize the lists dude!
      for(int S = 0;S < 2;++S)
         for(int i = 0;i < M/2;++i)
            for(int j = i + 1;j < M/2;++j)
               s2t[S][j][i] = s2t[S][i][j];

      //allocate
      _6j = new double * [2];

      for(int S = 0;S < 2;++S)
         _6j[S] = new double [2];

      //initialize
      _6j[0][0] = -0.5;
      _6j[0][1] = 0.5;
      _6j[1][0] = 0.5;
      _6j[1][1] = 1.0/6.0;

   }

   ++counter;

}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpm_c
 * if counter == 0, the lists containing the relationship between sp and tp basis.
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(TPM &tpm_c) : BlockMatrix(tpm_c){

   this->N = tpm_c.gN();
   this->M = tpm_c.gM();

   if(counter == 0){

      //allocatie van s2t
      s2t = new int ** [2];

      for(int i = 0;i < 2;++i){

         s2t[i] = new int * [M/2];
         s2t[i][0] = new int [M*M/4];

         for(int j = 1;j < M/2;++j)
            s2t[i][j] = s2t[i][j - 1] + M/2;

      }

      //allocatie van t2s
      t2s = new int ** [2];

      for(int i = 0;i < 2;++i){

         t2s[i] = new int * [this->gdim(i)];

         for(int j = 0;j < this->gdim(i);++j)
            t2s[i][j] = new int [2];

      }

      //initialisatie van de arrays
      int teller = 0;

      //symmetrical array: S = 0
      for(int a = 0;a < M/2;++a)
         for(int b = a;b < M/2;++b){

            s2t[0][a][b] = teller;

            t2s[0][teller][0] = a;
            t2s[0][teller][1] = b;

            ++teller;

         }

      //watch it!
      teller = 0;

      //antisymmetrical array: S = 1
      for(int a = 0;a < M/2;++a)
         for(int b = a + 1;b < M/2;++b){

            s2t[1][a][b] = teller;

            t2s[1][teller][0] = a;
            t2s[1][teller][1] = b;

            ++teller;

         }

      //symmetrize the lists dude!
      for(int S = 0;S < 2;++S)
         for(int i = 0;i < M/2;++i)
            for(int j = i + 1;j < M/2;++j)
               s2t[S][j][i] = s2t[S][i][j];

      //allocate
      _6j = new double * [2];

      for(int S = 0;S < 2;++S)
         _6j[S] = new double [2];

      //initialize
      _6j[0][0] = -0.5;
      _6j[0][1] = 0.5;
      _6j[1][0] = 0.5;
      _6j[1][1] = 1.0/6.0;

   }

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists t2s en s2t will be deleted.
 * 
 */
TPM::~TPM(){

   if(counter == 1){

      for(int S = 0;S < 2;++S){

         delete [] _6j[S];

         delete [] s2t[S][0];
         delete [] s2t[S];

         for(int i = 0;i < this->gdim(S);++i)
            delete [] t2s[S][i];

         delete [] t2s[S];

      }

      delete [] s2t;
      delete [] t2s;

      delete [] _6j;

   }

   --counter;

}

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The spinquantumnumber that identifies the block
 * @param a first sp index that forms the tp row index i of spin S, together with b
 * @param b second sp index that forms the tp row index i of spin S, together with a
 * @param c first sp index that forms the tp column index j of spin S, together with d
 * @param d second sp index that forms the tp column index j of spin S, together with c
 * @return the number on place TPM(i,j) with the right phase.
 */
double TPM::operator()(int S,int a,int b,int c,int d) const{

   if(S == 0){

      int i = s2t[0][a][b];
      int j = s2t[0][c][d];

      return (*this)(S,i,j);

   }
   else{

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[1][a][b];
         int j = s2t[1][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase*(*this)(S,i,j);

      }

   }

}

ostream &operator<<(ostream &output,TPM &tpm_p){

   for(int S = 0;S < 2;++S){

      output << S << "\t" << tpm_p.gdim(S) << "\t" << tpm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(S);++i)
         for(int j = 0;j < tpm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[S][i][0] << "\t" << tpm_p.t2s[S][i][1]

               << "\t" << tpm_p.t2s[S][j][0] << "\t" << tpm_p.t2s[S][j][1] << "\t" << tpm_p(S,i,j) << endl;

         }

      std::cout << std::endl;

   }

   return output;

}

/**
 * @return number of particles
 */
int TPM::gN(){

   return N;

}

/**
 * @return number of sp orbitals
 */
int TPM::gM(){

   return M;

}

/**
 * construct the spinsymmetrical hubbard hamiltonian with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int a,b,c,d;//sp (lattice sites here) orbitals

   double ward = 1.0/(N - 1.0);

   int sign;

   for(int S = 0;S < 2;++S){

      sign = 1 - 2*S;

      for(int i = 0;i < this->gdim(S);++i){

         a = t2s[S][i][0];
         b = t2s[S][i][1];

         for(int j = i;j < this->gdim(S);++j){

            c = t2s[S][j][0];
            d = t2s[S][j][1];

            (*this)(S,i,j) = 0;

            //eerst hopping
            if( (a == c) && ( ( (b + 1)%(M/2) == d ) || ( b == (d + 1)%(M/2) ) ) )
               (*this)(S,i,j) -= ward;

            if( (b == c) && ( ( (a + 1)%(M/2) == d ) || ( a == (d + 1)%(M/2) ) ) )
               (*this)(S,i,j) -= sign*ward;

            if( (a == d) && ( ( (b + 1)%(M/2) == c ) || ( b == (c + 1)%(M/2) ) ) )
               (*this)(S,i,j) -= sign*ward;

            if( (b == d) && ( ( (a + 1)%(M/2) == c ) || ( a == (c + 1)%(M/2) ) ) )
               (*this)(S,i,j) -= ward;

            //only on-site interaction for singlet tp states:
            if(S == 0)
               if(i == j && a == b)
                  (*this)(S,i,j) += 2.0*U;

            if(a == b)
               (*this)(S,i,j) /= std::sqrt(2.0);

            if(c == d)
               (*this)(S,i,j) /= std::sqrt(2.0);

         }
      }

   }

   this->symmetrize();

}

/**
 * The spincoupled Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,TPM &tpm_d){

   double a = 1;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * The spincoupled Q-like map: see primal-dual.pdf for more info (form: Q^S(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,TPM &tpm_d){

   //for inverse
   if(option == -1){

      B = (B*A + B*C*M - 2.0*C*C)/( A * (C*(M - 2.0) -  A) * ( A + B*M*(M - 1.0) - 2.0*C*(M - 1.0) ) );
      C = C/(A*(C*(M - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm(C,tpm_d);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   int sign;

   double norm;

   //loop over the spinblocks
   for(int S = 0;S < 2;++S){

      //symmetry or antisymmetry?
      sign = 1 - 2*S;

      for(int i = 0;i < this->gdim(S);++i){

         int a = t2s[S][i][0];
         int b = t2s[S][i][1];

         for(int j = i;j < this->gdim(S);++j){

            int c = t2s[S][j][0];
            int d = t2s[S][j][1];

            //determine the norm for the basisset
            norm = 1.0;

            if(S == 0){

               if(a == b)
                  norm /= std::sqrt(2.0);

               if(c == d)
                  norm /= std::sqrt(2.0);

            }

            //here starts the Q-map

            //the tp part
            (*this)(S,i,j) = A*tpm_d(S,i,j);

            //the np part
            if(i == j)
               (*this)(S,i,i) += ward;

            //and four sp parts:
            if(a == c)
               (*this)(S,i,j) -= norm*spm(b,d);

            if(b == c)
               (*this)(S,i,j) -= sign*norm*spm(a,d);

            if(a == d)
               (*this)(S,i,j) -= sign*norm*spm(b,c);

            if(b == d)
               (*this)(S,i,j) -= norm*spm(a,c);

         }
      }

   }

   this->symmetrize();

}

/**
 * initialize onto the TPM with S = 0 block empty and unit times trace in the S = 1 block
 */
void TPM::init(){

   double ward = N*(N - 1.0)/(M/2*(M/2 - 1.0))/3.0;


   for(int i = 0;i < gdim(0);++i)
      for(int j = 0;j < gdim(0);++j)
         (*this)(0,i,j) = 0.0;


   for(int i = 0;i < gdim(1);++i){

      (*this)(1,i,i) = ward;

      for(int j = i + 1;j < gdim(1);++j)
         (*this)(1,i,j) = (*this)(1,j,i) = 0.0;

   }

}

/**
 * orthogonal projection onto the space of traceless matrices of S_max (so with S = 0 with zero entries)
 */
void TPM::proj_E(){

   //first put the S = 0 block to zero
   (*this)[0] = 0;

   //then make the S = 1 block traceless
   double ward = 2 * this->trace()/(M/2*(M/2 - 1))/3.0;

   for(int i = 0;i < gdim(1);++i)
      (*this)(1,i,i) -= ward;

}

/**
 * Primal hessian map:\n\n
 * Hb = D_1 b D_1 + D_2 Q(b) D_2 + D_3 G(b) D_3 + D_4 T1(b) D_4 + D_5 T2(b) D5 \n\n
 * with D_1, D_2, D_3, D_4 and D_5 the P, Q, G, T1 and T2 blocks of the SUP D. 
 * @param b TPM domain matrix, hessian will act on it and the image will be put in this
 * @param D SUP matrix that defines the structure of the hessian map. (see primal-dual.pdf for more info)
 */
void TPM::H(TPM &b,SUP &D){

   this->L_map(D.tpm(0),b);

#ifdef __Q_CON

   //maak Q(b)
   TPM Qb(M,N);
   Qb.Q(1,b);

   TPM hulp(M,N);

   hulp.L_map(D.tpm(1),Qb);

   Qb.Q(1,hulp);

   *this += Qb;

#endif

#ifdef __G_CON

   //maak G(b)
   PHM Gb(M,N);
   Gb.G(b);

   PHM hulpje(M,N);

   hulpje.L_map(D.phm(),Gb);

   hulp.G(hulpje);

   *this += hulp;

#endif

#ifdef __T1_CON

   DPM T1b(M,N);
   T1b.T(b);

   DPM hulp_T1(M,N);

   hulp_T1.L_map(D.dpm(),T1b);

   hulp.T(hulp_T1);

   *this += hulp;

#endif

#ifdef __T2_CON

   PPHM T2b(M,N);
   T2b.T(b);

   PPHM hulp_T2(M,N);

   hulp_T2.L_map(D.pphm(),T2b);

   hulp.T(hulp_T2);

   *this += hulp;

#endif

   this->proj_E();

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */
void TPM::S(int option,TPM &tpm_d){

   double a = 1.0;
   double b = 0.0;
   double c = 0.0;

#ifdef __Q_CON

   a += 1.0;
   b += (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   c += (2.0*N - M)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __G_CON

   a += 4.0;
   c += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));

#endif

#ifdef __T1_CON

   a += M - 4.0;
   b += (M*M*M - 6.0*M*M*N -3.0*M*M + 12.0*M*N*N + 12.0*M*N + 2.0*M - 18.0*N*N - 6.0*N*N*N)/( 3.0*N*N*(N - 1.0)*(N - 1.0) );
   c -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );

#endif

#ifdef __T2_CON

   a += 5.0*M - 8.0;
   b += 2.0/(N - 1.0);
   c += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));

#endif

   this->Q(option,a,b,c,tpm_d);

}

/**
 * Collaps a SUP matrix S onto a TPM matrix like this:\n\n
 * sum_i Tr (S u^i)f^i = this
 * @param option = 0, project onto full symmetric matrix space, = 1 project onto traceless symmetric matrix space
 * @param S input SUP
 */
void TPM::collaps(int option,SUP &S){

   *this = S.tpm(0);

   TPM hulp(M,N);

   hulp.Q(1,S.tpm(1));

   *this += hulp;

#ifdef __G_CON

   hulp.G(S.phm());

   *this += hulp;

#endif

#ifdef __T1_CON

   hulp.T(S.dpm());

   *this += hulp;

#endif

#ifdef __T2_CON

   hulp.T(S.pphm());

   *this += hulp;

#endif

   if(option == 1)
      this->proj_E();

}

/** 
 * Construct the pairing hamiltonian with a single particle spectrum
 * @param pair_coupling The strenght of the pairing interaction
 */
void TPM::sp_pairing(double pair_coupling){

   double *E = new double [M/2];

   //single particle spectrum
   for(int a = -M/2;a < 0;++a)
      E[M/2 + a] = (double) a;

   double *x = new double [M/2];

   //pairing interaction term
   for(int a = 0;a < M/2;++a)
      x[a] = 1.0;

   //normeren op 1/2
   double ward = 0.0;

   for(int i = 0;i < M/2;++i)
      ward += x[i]*x[i];

   ward *= 2.0;

   for(int a = 0;a < M/2;++a)
      x[a] /= std::sqrt(ward);

   int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < this->gdim(S);++i){

         a = t2s[S][i][0];
         b = t2s[S][i][1];

         //sp stuk
         (*this)(S,i,i) = (E[a] + E[b])/(N - 1.0);

         if(a == b)
            (*this)(S,i,i) -= 2.0*pair_coupling*x[a]*x[a];

         for(int j = i + 1;j < this->gdim(S);++j){

            c = t2s[S][j][0];
            d = t2s[S][j][1];

            if(a == b && c == d)
               (*this)(S,i,j) = -2.0*pair_coupling*x[a]*x[c];

         }

      }

   }

   this->symmetrize();

   delete [] E;

}

/**
 * Will print out a spin uncoupled version of the TPM (*this) in uncoupled sp coordinates to the file with name, filename
 */
void TPM::uncouple(const char *filename){

   ofstream output(filename);

   output.precision(10);

   int a,b,c,d;
   int s_a,s_b,s_c,s_d;

   double ward_0;
   double ward_1;

   //the uncoupled row vector i
   for(int alpha = 0;alpha < M;++alpha)
      for(int beta = alpha + 1;beta < M;++beta){

         //watch it now, 0 is down, 1 is up.
         a = alpha/2;
         s_a = alpha%2;

         b = beta/2;
         s_b = beta%2;

         //the collom vector j
         for(int gamma = alpha;gamma < M;++gamma)
            for(int delta = gamma + 1;delta < M;++delta){

               c = gamma/2;
               s_c = gamma%2;

               d = delta/2;
               s_d = delta%2;

               //eerst S = 0 stuk
               if(s_a != s_b && s_c != s_d){

                  if(s_a == s_c)
                     ward_0 = 0.5 * (*this)(0,a,b,c,d);
                  else
                     ward_0 = -0.5 * (*this)(0,a,b,c,d);

                  //normering
                  if(a == b)
                     ward_0 *= std::sqrt(2.0);

                  if(c == d)
                     ward_0 *= std::sqrt(2.0);

               }
               else
                  ward_0 = 0.0;

               //dan S = 1 stuk
               if( (s_a == s_b) && (s_c == s_d) && (s_a == s_c) )
                  ward_1 = (*this)(1,a,b,c,d);
               else if(s_a != s_b && s_c != s_d)
                  ward_1 = 0.5*(*this)(1,a,b,c,d);
               else
                  ward_1 = 0.0;

               output << alpha << "\t" << beta << "\t" << gamma << "\t" << delta << "\t" << ward_0 + ward_1 << std::endl;

            }

      }

}

/**
 * The G down map, maps a PHM object onto a TPM object using the G map.
 * @param phm input PHM
 */
void TPM::G(PHM &phm){

   SPM spm(1.0/(N - 1.0),phm);

   int sign;
   int a,b,c,d;

   for(int S = 0;S < 2;++S){

      sign = 1 - 2*S;

      for(int i = 0;i < this->gdim(S);++i){

         a = t2s[S][i][0];
         b = t2s[S][i][1];

         for(int j = i;j < this->gdim(S);++j){

            c = t2s[S][j][0];
            d = t2s[S][j][1];

            //init
            (*this)(S,i,j) = 0.0;

            //ph part
            for(int Z = 0;Z < 2;++Z)
               (*this)(S,i,j) -= this->gdeg(Z)*_6j[S][Z] * ( phm(Z,a,d,c,b) + phm(Z,b,c,d,a) + sign*phm(Z,b,d,c,a) + sign*phm(Z,a,c,d,b) );

            //4 sp parts
            if(b == d)
               (*this)(S,i,j) += spm(a,c);

            if(a == c)
               (*this)(S,i,j) += spm(b,d);

            if(a == d)
               (*this)(S,i,j) += sign*spm(b,c);

            if(b == c)
               (*this)(S,i,j) += sign*spm(a,d);

            //norm of the basisset:
            if(a == b)
               (*this)(S,i,j) /= std::sqrt(2.0);

            if(c == d)
               (*this)(S,i,j) /= std::sqrt(2.0);

         }

      }

   }

   this->symmetrize();

}

/**
 * Construct a spincoupled TPM matrix out of a spincoupled DPM matrix, for the definition and derivation see symmetry.pdf
 * @param dpm input DPM
 */
void TPM::bar(DPM &dpm){

   int a,b,c,d;

   double ward;

   //first the S = 0 part, easiest:
   for(int i = 0;i < this->gdim(0);++i){

      a = t2s[0][i][0];
      b = t2s[0][i][1];

      for(int j = i;j < this->gdim(0);++j){

         c = t2s[0][j][0];
         d = t2s[0][j][1];

         (*this)(0,i,j) = 0.0;

         //only total S = 1/2 can remain because cannot couple to S = 3/2 with intermediate S = 0
         for(int l = 0;l < M/2;++l)
            (*this)(0,i,j) += 2.0 * dpm(0,0,a,b,l,0,c,d,l);

      }
   }

   //then the S = 1 part:
   for(int i = 0;i < this->gdim(1);++i){

      a = t2s[1][i][0];
      b = t2s[1][i][1];

      for(int j = i;j < this->gdim(1);++j){

         c = t2s[1][j][0];
         d = t2s[1][j][1];

         (*this)(1,i,j) = 0.0;

         for(int Z = 0;Z < 2;++Z){//loop over the dpm blocks: S = 1/2 and 3/2 = Z + 1/2

            ward = 0.0;

            for(int l = 0;l < M/2;++l)
               ward += dpm(Z,1,a,b,l,1,c,d,l);

            ward *= (2 * (Z + 0.5) + 1.0)/3.0;

            (*this)(1,i,j) += ward;

         }

      }
   }

   this->symmetrize();

}

/** 
 * The T1-down map that maps a DPM on TPM. This is just a Q-like map using the TPM::bar (dpm) as input.
 * @param dpm the input DPM matrix
 */
void TPM::T(DPM &dpm){

   TPM tpm(M,N);
   tpm.bar(dpm);

   double a = 1;
   double b = 1.0/(3.0*N*(N - 1.0));
   double c = 0.5/(N - 1.0);

   this->Q(1,a,b,c,tpm);

}

/**
 * The bar function that maps a PPHM object onto a TPM object by tracing away the last pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void TPM::bar(PPHM &pphm){

   int a,b,c,d;

   double ward;

   for(int Z = 0;Z < 2;++Z){

      for(int i = 0;i < this->gdim(Z);++i){

         a = t2s[Z][i][0];
         b = t2s[Z][i][1];

         for(int j = i;j < this->gdim(Z);++j){

            c = t2s[Z][j][0];
            d = t2s[Z][j][1];

            (*this)(Z,i,j) = 0.0;

            for(int S = 0;S < 2;++S){//loop over three particle spin: 1/2 and 3/2

               ward = (2.0*(S + 0.5) + 1.0)/(2.0*Z + 1.0);

               for(int l = 0;l < M/2;++l)
                  (*this)(Z,i,j) += ward * pphm(S,Z,a,b,l,Z,c,d,l);

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * The spincoupled T2-down map that maps a PPHM on a TPM object.
 * @param pphm Input PPHM object
 */
void TPM::T(PPHM &pphm){

   //first make the bar tpm
   TPM tpm(M,N);
   tpm.bar(pphm);

   //then make the bar phm
   PHM phm(M,N);
   phm.bar(pphm);

   //also make the bar spm with the correct scale factor
   SPM spm(M,N);
   spm.bar(0.5/(N - 1.0),pphm);

   int a,b,c,d;
   int sign;

   double norm;

   for(int S = 0;S < 2;++S){

      sign = 1 - 2*S;

      for(int i = 0;i < this->gdim(S);++i){

         a = t2s[S][i][0];
         b = t2s[S][i][1];

         for(int j = i;j < this->gdim(S);++j){

            c = t2s[S][j][0];
            d = t2s[S][j][1];

            //determine the norm for the basisset
            norm = 1.0;

            if(S == 0){

               if(a == b)
                  norm /= std::sqrt(2.0);

               if(c == d)
                  norm /= std::sqrt(2.0);

            }

            //first the tp part
            (*this)(S,i,j) = tpm(S,i,j);

            //sp part, 4 terms:
            if(b == d)
               (*this)(S,i,j) += norm * spm(a,c);

            if(a == d)
               (*this)(S,i,j) += sign * norm * spm(b,c);

            if(b == c)
               (*this)(S,i,j) += sign * norm * spm(a,d);

            if(a == c)
               (*this)(S,i,j) += norm * spm(b,d);

            //ph part:
            for(int Z = 0;Z < 2;++Z)
               (*this)(S,i,j) -= norm * (2.0 * Z + 1.0) * _6j[S][Z] * ( phm(Z,d,a,b,c) + sign * phm(Z,d,b,a,c) + sign * phm(Z,c,a,b,d) + phm(Z,c,b,a,d) );

         }
      }

   }

   this->symmetrize();

}

/**
 * @return The expectation value of the total spin for the TPM.
 */
double TPM::spin(){

   double ward = 0.0;

   for(int i = 0;i < this->gdim(0);++i)
      ward += -1.5 * (N - 2.0)/(N - 1.0) * (*this)(0,i,i);

   for(int i = 0;i < this->gdim(1);++i)
      ward += 3.0 * ( -1.5 * (N - 2.0)/(N - 1.0) + 2.0 ) * (*this)(1,i,i);

   return ward;

}

/**
 * Construct the right hand side of the Newton equation for the determination of the search direction, 
 * the gradient of the potential:
 * @param t scaling factor of the potential
 * @param ham Hamiltonian of the current problem
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 */
void TPM::constr_grad(double t,TPM &ham,SUP &P){

   //eerst P conditie 
   *this = P.tpm(0);

   //de Q conditie toevoegen

#ifdef __Q_CON

   TPM hulp(M,N);

   hulp.Q(1,P.tpm(1));

   *this += hulp;

#endif

   //de G conditie indien nodig

#ifdef __G_CON

   hulp.G(P.phm());

   *this += hulp;

#endif

   //de T_1 conditie toevoegen indien nodig

#ifdef __T1_CON

   hulp.T(P.dpm());

   *this += hulp;

#endif

#ifdef __T2_CON

   hulp.T(P.pphm());

   *this +=hulp;
#endif

   this->dscal(t);

   *this -= ham;

   this->proj_E();

}

/**
 * solve the Newton equations for the determination of the search direction,
 * @param t scaling factor of the potential
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 * @param b right hand side (the gradient constructed int TPM::constr_grad)
 * @return nr of iterations needed to converge to the desired accuracy
 */
int TPM::solve(double t,SUP &P,TPM &b){

   int iter = 0;

   //delta = 0
   *this = 0;

   //residu:
   TPM r(b);

   //norm van het residu
   double rr = r.ddot(r);

   //enkele variabelen
   double rr_old,ward;

   TPM Hb(M,N);

   while(rr > 1.0e-7){ 

      Hb.H(t,b,P);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe variabelen berekenen en oude overdragen
      rr_old = rr;
      rr = r.ddot(r);

      //nieuwe b nog:
      b.dscal(rr/rr_old);

      b += r;

      ++iter;

   }

   return iter;

}

/**
 * perform a line search what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,SUP &P,TPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta(M,N);

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp(M,N);

   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * The hessian-map of the Newton system:
 * @param t potential scaling factor
 * @param b the TPM on which the hamiltonian will work, the image will be put in (*this)
 * @param P the SUP matrix containing the constraints, (can be seen as the metric).
 */
void TPM::H(double t,TPM &b,SUP &P){

   //eerst de P conditie:

   this->L_map(P.tpm(0),b);

   //de Q conditie toevoegen:
#ifdef __Q_CON

   //hulpje
   TPM hulp(M,N);

   //maak Q(b)
   TPM Q_b(M,N);
   Q_b.Q(1,b);

   //stop Q(rdm)^{-1}Q(b)Q(rdm)^{-1} in hulp
   hulp.L_map(P.tpm(1),Q_b);

   //maak Q(hulp) en stop in Q_b
   Q_b.Q(1,hulp);

   //en tel op bij this
   *this += Q_b;

#endif

   //de G conditie toevoegen:

#ifdef __G_CON

   //hulpje voor het PHM stuk
   PHM hulp_ph(M,N);
   PHM G_b(M,N);

   //stop G(b) in G_b
   G_b.G(b);

   //bereken G(rdm)^{-1}G(b)G(rdm)^{-1} en stop in hulp_ph
   hulp_ph.L_map(P.phm(),G_b);

   //tenslotte nog de antisymmetrische G hierop:
   hulp.G(hulp_ph);

   //en optellen bij this
   *this += hulp;

#endif

   //de T_1 conditie toevoegen

#ifdef __T1_CON

   //hulpjes voor het DPM stuk
   DPM hulp_dp(M,N);
   DPM T1_b(M,N);

   //stop T1(b) in T1_b
   T1_b.T(b);

   hulp_dp.L_map(P.dpm(),T1_b);

   hulp.T(hulp_dp);

   *this += hulp;

#endif

#ifdef __T2_CON

   PPHM hulp_pph(M,N);
   PPHM T2_b(M,N);

   T2_b.T(b);

   hulp_pph.L_map(P.pphm(),T2_b);

   hulp.T(hulp_pph);

   *this+=hulp;

#endif

   //nog schalen met t:
   this->dscal(t);

   //en projecteren op spoorloze ruimte
   this->proj_E();

}

/**
 * perform a line search what step size in along the Newton direction is ideal, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param rdm TPM containing the current approximation of the rdm.
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,TPM &rdm,TPM &ham){

   SUP P(M,N);

   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham);

}

/**
 * method specific for the S = MAX program, inverts only the S = 1 block
 */
void TPM::pseudo_invert(){

   (*this)[1].invert();

}

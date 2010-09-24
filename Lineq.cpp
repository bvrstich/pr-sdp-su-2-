#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * standard constructor, only norm constraint:
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
Lineq::Lineq(int M,int N){

   this->nr = 1;

   this->M = M;
   this->N = N;

   //allocate the bitch.
   allocate();

   //initialize the constraint
   E[0]->set_unit();

   e[0] = N*(N - 1.0)/2.0;

   orthogonalize();//speaks for itself, doesn't it?

}

/**
 * constructor, norm and S^2 constraint:
 * @param M nr of sp orbitals
 * @param N nr of particles
 * @param spin the size of the spin you want to project on: e[1]
 */
Lineq::Lineq(int M,int N,double spin){

   this->nr = 2;

   this->M = M;
   this->N = N;

   //allocate the bitch.
   allocate();

   //initialize the constraints:first the norm
   E[0]->set_unit();

   e[0] = N*(N - 1.0)/2.0;

   //then the spin in the "1" block
   E[1]->set_S_2();

   e[1] = spin * ( spin + 1.0 );

   orthogonalize();//speaks for itself, doesn't it?

}

/**
 * copy constructor:
 * @param lineq the Lineq object that will be copied into this
 */
Lineq::Lineq(const Lineq &lineq){

   this->nr = lineq.gnr();

   this->M = lineq.gM();
   this->N = lineq.gN();

   E = new TPM * [nr];
   e = new double [nr];

   E_ortho = new TPM * [nr];
   e_ortho = new double [nr];

   //copy the that which needs to be copied
   for(int i = 0;i < nr;++i){

      E[i] = new TPM(lineq.gE(i)); 
      e[i] = lineq.ge(i);

      //no need to orthogonalize, just copy!
      E_ortho[i] = new TPM(lineq.gE_ortho(i)); 
      e_ortho[i] = lineq.ge_ortho(i);

   }

}

/**
 * destructor
 */
Lineq::~Lineq(){

   for(int i = 0;i < nr;++i){

      delete E[i];
      delete E_ortho[i];

   }

   delete [] E;
   delete [] E_ortho;

   delete [] e;
   delete [] e_ortho;

}

/**
 * @return nr of particles
 */
int Lineq::gN() const {

   return N;

}

/**
 * @return nr of constraints
 */
int Lineq::gnr() const {

   return nr;

}

/**
 * @return nr of sp orbitals
 */
int Lineq::gM() const {

   return M;

}

/**
 * access to the individual constraint TPM's
 * @param i the index
 * @return The E TPM on index i.
 */
TPM &Lineq::gE(int i) const {

   return *E[i];

}

/**
 * access to the individual constraint values
 * @param i the index
 * @return the e values on index i: e[i] or something
 */
double &Lineq::ge(int i) const{

   return e[i];

}

/**
 * access to the individual orthogonalized constraint TPM's: private function
 * @param i the index
 * @return The E_ortho TPM on index i.
 */
TPM &Lineq::gE_ortho(int i) const {

   return *E_ortho[i];

}

/**
 * access to the individual orthogonalized constraint values: private function
 * @param i the index
 * @return the e values on index i: e_ortho[i] or something
 */
double &Lineq::ge_ortho(int i) const{

   return e_ortho[i];

}

ostream &operator<<(ostream &output,Lineq &lineq_p){

   output << "first print the constraint matrices:";
   output << endl;
   output << endl;

   for(int i = 0;i < lineq_p.gnr();++i){

      output << "constraint nr :" << i << endl;
      output << endl;

      output << lineq_p.gE(i);

   }

   output << endl;
   output << endl;
   output << "the desired values are:" << endl;
   output << endl;

   for(int i = 0;i < lineq_p.gnr();++i)
      output << i << "\t" << lineq_p.ge(i) << endl;

   return output;

}

/**
 * orthogonalize the constraints, will take E and e, and construct E_ortho and e_ortho with them.
 */
void Lineq::orthogonalize(){

   //construct the overlapmatrix of the E's
   Matrix S(nr);

   for(int i = 0;i < nr;++i)
      for(int j = i;j < nr;++j)
         S(i,j) = E[i]->ddot(*E[j]);

   S.symmetrize();

   //take the inverse square root
   S.sqrt(-1);

   //make the orthogonal ones:
   for(int i = 0;i < nr;++i){

      *E_ortho[i] = 0;
      e_ortho[i] = 0;

      for(int j = 0;j < nr;++j){

         E_ortho[i]->daxpy(S(i,j),*E[j]);
         e_ortho[i] += S(i,j) * e[j];

      }

   }

}

/**
 * allocate the array's and memory needed for the program, reusable code for the different constructors
 * private because I say it's private.
 */
void Lineq::allocate(){

   //the regular ones
   E = new TPM * [nr];
   e = new double [nr];

   for(int i = 0;i < nr;++i)
      E[i] = new TPM(M,N);

   //the orthogonal ones
   E_ortho = new TPM * [nr];
   e_ortho = new double [nr];

   for(int i = 0;i < nr;++i)
      E_ortho[i] = new TPM(M,N);

}

/**
 * @mainpage 
 * This is an implementation of the spinsymmetrical, dual only, potential reduction interior point method
 * for optimizing the second order density matrix using the P, Q, G, T1 and T2 N-representability conditions.
 * Compiling can be done with the options PQ, PQG, PQGT1, PQGT2 and PQGT (for all conditions active) with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 22-02-2010
 */

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity density matrix normed on the particle number and minimize the 
 * ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void){

   cout.precision(10);

   const int M = 8;//dim sp hilbert space
   const int N = 4;//nr of particles

   //linear equalities object
   Lineq lineq(M,N,1.0);

   //hamiltoniaan
   TPM ham(M,N);

   ham.hubbard(1.0);

   TPM rdm(M,N);
   rdm.init(lineq);

   TPM backup_rdm(rdm);

   double t = 1.0;
   double tolerance = 1.0e-5;

   //outer iteration: scaling of the potential barrier
   while(t > 1.0e-12){

      cout << t << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) << endl;

      double convergence = 1.0;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      while(convergence > tolerance){

         SUP P(M,N);

         P.fill(rdm);

         P.invert();

         //eerst -gradient aanmaken:
         TPM grad(M,N);

         grad.constr_grad(t,ham,P,lineq);

         //dit wordt de stap:
         TPM delta(M,N);

         //los het hessiaan stelsel op:
         cout << delta.solve(t,P,grad,lineq) << endl;

         //line search
         double a = delta.line_search(t,P,ham);

         //rdm += a*delta;
         rdm.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

      }

      t /= 2.0;

      //what is the tolerance for the newton method?
      tolerance = 1.0e-5*t;

      if(tolerance < 1.0e-12)
         tolerance = 1.0e-12;

      //extrapolatie:
      TPM extrapol(rdm);

      extrapol -= backup_rdm;

      //overzetten voor volgende stap
      backup_rdm = rdm;

      double a = extrapol.line_search(t,rdm,ham);

      rdm.daxpy(a,extrapol);

   }

   cout << endl;
   cout << "Energy = \t" << rdm.ddot(ham) << endl;
   cout << endl;
   cout << "S^2 = \t" << rdm.spin() << endl;

   return 0;
}

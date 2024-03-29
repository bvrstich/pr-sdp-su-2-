#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

class SUP;
class PHM;
class DPM;
class PPHM;

/**
 * @author Brecht Verstichel
 * @date 19-04-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry included, it inherits alle the function from its mother 
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,TPM &tpm_p);

   public:
      
      //constructor
      TPM(int M,int N);

      //copy constructor
      TPM(TPM &);

      //destructor
      virtual ~TPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode and with spin quantumnumer
      double operator()(int S,int a,int b,int c,int d) const;

      //geef N terug
      int gN();

      //geef M terug
      int gM();

      void hubbard(double U);

      //Q afbeelding en zijn inverse
      void Q(int option,TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,TPM &);

      void unit();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(TPM &b,SUP &D);

      //los het stelsel op
      int solve(TPM &b,SUP &D);

      void min_unit(double scale);

      void min_qunit(double scale);

      void collaps(int option,SUP &);

      void sp_pairing(double );

      void uncouple(const char *);

      //G down afbeelding
      void G(PHM &);

      //trace one pair of indices of DPM
      void bar(DPM &);

      //T1 down
      void T(DPM &);

      //trace last pair of indices of PPHM
      void bar(PPHM &);

      //T2 down
      void T(PPHM &);

      //return the spin
      double spin();

      void constr_grad(double t,TPM &,SUP &);

      int solve(double t,SUP &,TPM &);

      double line_search(double t,SUP &P,TPM &ham);

      double line_search(double t,TPM &,TPM &);

      void H(double t,TPM &b,SUP &P);

      void constr_sp_diag(int);


   private:

      //!static list of dimension [2][dim[i]][2] that takes in a tp index i and a spinquantumnumber S, and returns two sp indices: a = t2s[S][i][0] and b = t2s[S][i][1]
      static int ***t2s;

      //!static list of dimension [M][M] that takes two sp indices a,b and a spinquantumnumber S, and returns a tp index i: i = s2t[S][a][b]
      static int ***s2t;

      //!list of 6j symbols needed.
      static double **_6j;

      //!static counter that counts the number of TPM objects running in the program
      static int counter;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

};

#endif

#ifndef LINEQ_H
#define LINEQ_H

#include <iostream>
#include <cstdlib>

using std::ostream;

#include "TPM.h"
#include "SUP.h"

/**
 * @author Brecht Verstichel
 * @date 09-09-2010\n\n
 * This is a class written for handling linear equality constraints in the program. Contains a collection of the matrices
 * on which the tpm has to have a certain projection, these values are stored in a pointer of doubles.
 * 
 */

class Lineq {

   /**
    * Output stream operator overloaded: prints all the constraint matrices with a whitespace between them.
    * @param output The stream to which you are writing (e.g. cout)
    * @param lineq_p The Lineq you want to print
    */
   friend ostream &operator<<(ostream &output,Lineq &lineq_p);

   public:

      //norm constructor
      Lineq(int M,int N);

      //norm + spinsize constructor
      Lineq(int M,int N,double spin);

      //copy constructor
      Lineq(const Lineq &);

      //destructor
      virtual ~Lineq();

      int gN() const;

      int gnr() const;

      int gM() const;

      TPM &gE(int) const;

      double &ge(int) const;

      TPM &gE_ortho(int) const;

      double &ge_ortho(int) const;

   private:

      void allocate();

      void orthogonalize();

      //!double pointer to TPM object, will contain the linear equality constraints
      TPM **E;

      //!pointer of doubles, will contain the values of the projections. (the desired equalities)
      double *e;

      //!orthogonalized constraints, these will be hidden from the public.
      TPM **E_ortho;

      //!the values accompanying the orthogonalized constraints
      double *e_ortho;

      //!nr of contstraints
      int nr;

      //!nr of particles
      int N;

      //!nr of sp orbs
      int M;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/

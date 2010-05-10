#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <iostream>
#include <cstdlib>

#include "Matrix.h"

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 15-04-2010\n\n
 * This is a class written for symmetric block matrices. It containts an array of Matrix objects and an 
 * array containing the dimensions of the different block. It redefines all the member functions of the Matrix
 * class, which uses the lapack and blas routines for matrix computations.
 */

class BlockMatrix{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << blockmatrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << blockmatrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param blockmatrix_p the BlockMatrix you want to print
    */
   friend ostream &operator<<(ostream &output,BlockMatrix &blockmatrix_p);

   public:

      //constructor
      BlockMatrix(int n);

      //copy constructor
      BlockMatrix(BlockMatrix &);

      //destructor
      virtual ~BlockMatrix();

      void setMatrixDim(int,int,int);

      Matrix &operator[](int);

      //overload equality operator
      BlockMatrix &operator=(BlockMatrix &);

      BlockMatrix &operator=(double );

      //overload += operator
      BlockMatrix &operator+=(BlockMatrix &);

      //overload -= operator
      BlockMatrix &operator-=(BlockMatrix &);

      BlockMatrix &daxpy(double alpha,BlockMatrix &);

      BlockMatrix &operator/=(double );

      BlockMatrix &mprod(BlockMatrix &,BlockMatrix &);

      //easy to change the numbers
      double &operator()(int block,int i,int j);

      //easy to access the numbers
      double operator()(int block,int i,int j) const;

      int gnr();

      int gdim(int);

      int gdeg(int);

      double trace();

      //void diagonalize(double **eigenvalues);

      double ddot(BlockMatrix &);

      void invert();

      void dscal(double alpha);

      void fill_Random();

      //positieve of negatieve vierkantswortel uit de matrix
      void sqrt(int option);

      void mdiag(double *diag);

      void L_map(BlockMatrix &,BlockMatrix &);

      void symmetrize();

      void out(const char *);

   private:

      //!pointer to Matrix objects, will contain the different blocks
      Matrix **blockmatrix;

      //!nr of blocks in the blockmatrix
      int nr;

      //!array of int's containing the dimensions of the different blocks
      int *dim;

      //!array of flags containing the information if the memory for the blockmatrices has been allocated (flag = 1) or not (flag = 0)
      int *flag;

      //!degeneracy of the blocks
      int *degen;

};

#endif

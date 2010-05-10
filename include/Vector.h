#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;

#include "lapack.h"
#include "Matrix.h"

class Vector;

/**
 * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
 * ifstream object and type:\n\n
 * object << vector << endl;\n\n
 * For output onto the screen type: \n\n
 * cout << vector << endl;\n\n
 * @param output The stream to which you are writing (e.g. cout)
 * @param vector_p de Vector you want to print
 */
ostream &operator<<(ostream &output,Vector &vector_p);

/**
 * @author Brecht Verstichel
 * @date 15-04-2010\n\n
 * This is a class written for vectors. It will contain the eigenvalues of the TPM, etc. Matrices. It is a template class,
 * corresponding to the different VectorType's that can be put in, it will automatically get the right dimension.
 * It is a wrapper around a pointer and redefines much used lapack and blas routines as memberfunctions
 */
class Vector{

   public:

      //construct with as input a dimension
      Vector(int );

      //construct with as input a Matrix
      Vector(Matrix &);

      //copy constructor
      Vector(Vector &);

      //destructor
      virtual ~Vector();

      //overload equality operator
      Vector &operator=(Vector &);

      Vector &operator=(double );

      //overload += operator
      Vector &operator+=(Vector &);

      //overload -= operator
      Vector &operator-=(Vector &);

      Vector &daxpy(double alpha,Vector &);

      Vector &operator/=(double );

      void diagonalize(Matrix &);

      //easy to change the numbers
      double &operator[](int i);

      //easy to access the numbers
      double operator[](int i) const;

      void fill_Random();

      //get the pointer to the vector
      double *gVector();

      int gn();

      double sum();

      double log_product();

      double ddot(Vector &);

      void dscal(double alpha);

      double min();
      
      double max();

      double centerpot(double );

      double lsfunc(double);

   private:

      //!pointer of doubles, contains the numbers, the vector
      double *vector;

      //!dimension of the vector
      int n;

};

#endif

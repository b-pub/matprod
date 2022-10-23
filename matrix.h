/*
 * Operand-collecting Matrix multiplication.
 *
 * Copyright (c) 2016 Brent Burton
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
 * KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 * WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
/*
 * MatrixProduct calculates a chained matrix multiplication.
 *
 * The way this works is based on the way C++ parses
 * an expression like "a * b * c" where a,b,c are Matrix
 * instances. The first multiplication operator*(Matrix,Matrix)
 * creates and returns a MatrixProduct with the first two
 * matrices collected (a and b). The next multiplication is
 * via operator*(MatrixProduct,Matrix) which simply appends
 * Matrix c into the collection. Any other multiplied matrices
 * are collected in the same way.
 *
 * The MatrixProduct does not know when the chain ends. Since
 * expressions look like "Matrix result = a * b * c;" the final
 * step is the Matrix::operator=(MatrixProduct). This is a sure
 * point that the chain is complete, and the entire chain can be
 * calculated. Thus, M::operator=(MP) calls
 * MatrixProduct::calculate() which plans and performs the
 * calculations.
 *
 * The goal of this creation was to explore a way to do delayed
 * computation. This allows arbitrary optimizations to be
 * performed inside the ::calculate() method. The current
 * algorithm is the Matrix Chain Multiplication
 * dynamic-programming example from
 * "Introduction to Algorithms", by Cormen, Leiserson, and Rivest, 1994.
 */
#ifndef _MATRIX_H
#define _MATRIX_H

#include <string>
#include <vector>
#include <limits>

namespace bb {

    class Matrix;

/*
 * MatrixProduct calculates the product of a chain of matrices.
 * This class isn't intended to be used directly. The various
 * operators defined elsewhere manipulate MatrixProducts.
 */
    class MatrixProduct
    {
        std::vector<Matrix> m_matrices;
        size_t m_numMults;

      public:

        /**
         * Create a MatrixProduct with initial two matrices.
         */
        MatrixProduct(const Matrix &m1, const Matrix &m2);

        /**
         * Accumulate another matrix for the product.
         */
        MatrixProduct& operator*(const Matrix &m);

        /**
         * Perform calculation planning and computation. In this
         * version, a boolean value is returned to report whether the
         * chain is compatible w.r.t. matrix dimensions.
         *   If this was implemented with exceptions, the return type
         * here would be Matrix, and if the chain was not compatible
         * an appropriate exception should be thrown.
         */
        bool calculate(Matrix &OUTresult);

        size_t getNaiveCost(); // Returns #mults for naive serial computation
        size_t getOptimizedCost(); // Returns #mults for optimized calculation

      private:

        void makeChainOrder(Matrix &OUTresult);
        Matrix matrixChainMultiply(int **s, int i, int j);
    };

/**
 * A 2D Matrix of floats
 */
    class Matrix
    {
      private:
        /**
         * MatrixPriv holds private matrix internals, possibly shared
         * across many Matrix objects via reference counting.
         */
        struct MatrixPriv {
          private:
            uint refs;                  // ref counted shared structure.
          public:
            std::string name;
            int rows, cols;
            float **mat;

            MatrixPriv();
            MatrixPriv(int m, int n, std::string _name);

            MatrixPriv* addRef();
            void release();

          private:
            ~MatrixPriv();
        };
        MatrixPriv *m_priv;
        bool m_ok;
        void init(int m, int n, std::string name);

        // Store costs of each multiplication style for comparison:
        size_t m_naiveCost;
        size_t m_optimizedCost;

      public:

        Matrix();                                 // new empty matrix
        Matrix(int m, int n, std::string name); // alloc a new
        Matrix(const MatrixProduct& mp);   // new matrix from product
        Matrix(const Matrix& m);

        bool ok() const;                        // is Matrix state OK?

        std::string name() const;
        void setName(std::string);

        int rows() const;
        int cols() const;

        size_t getNaiveCost() const     { return m_naiveCost; }
        size_t getOptimizedCost() const { return m_optimizedCost; }

        /**
         * Matrix data accessor. Access elements as M[i][j],
         * and i,j are not fully range-checked.
         */
        float* operator[](int);

        /**
         * Replace *this' guts with those shared from o.
         *
         * NOTE that this may not be the correct semantic, and in
         * spite of the const, the private data can still be
         * modified by the receiver.
         */
        Matrix& operator=(const Matrix& o);

        /**
         * Multiply *this against m, returning a new Matrix with the
         * result. The numMults parameter is set to the number of
         * multiplications required to calculate the product.
         */
        Matrix mult(Matrix& m, size_t& numMults);

        ~Matrix();
    };

}

/*
 * Global operator support to kick things off.
 */
bb::MatrixProduct operator*(const bb::Matrix& m1,
                            const bb::Matrix& m2);

#endif // _MATRIX_H

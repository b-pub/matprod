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
#include <iostream>
#include <iomanip>
#include <cassert>
#include "matrix.h"

// ---------------- MatrixPriv

bb::Matrix::MatrixPriv::MatrixPriv()
    : refs(1)
    , name("empty")
    , rows(0), cols(0), mat(0)
{ }

bb::Matrix::MatrixPriv::MatrixPriv(int m, int n, std::string _name)
    : refs(1)
    , name(_name)
    , rows(0), cols(0), mat(0)
{
    if (m>0 && n>0) {
        rows = m;
        cols = n;
        mat = new float*[rows];
        for (uint r=0; r<rows; r++)
            mat[r] = new float[cols];
    }
}

bb::Matrix::MatrixPriv::~MatrixPriv()
{
    if (rows && cols)
    {
        for (uint r=0; r<rows; r++)
            delete[] mat[r];
        delete[] mat;
    }
}

bb::Matrix::MatrixPriv* bb::Matrix::MatrixPriv::addRef()
{
    refs++;
    return this;
}

void bb::Matrix::MatrixPriv::release()
{
    assert(refs > 0);
    if (--refs == 0)
        delete this;
}

// ---------------- Matrix

bb::Matrix::Matrix(const MatrixProduct& mp)
    : m_priv(0)
    , m_ok(false)
{
    Matrix result;
    MatrixProduct &ncmp = const_cast<MatrixProduct&>(mp);
    m_ok = ncmp.calculate(result);

    if (m_ok) {
        m_priv = result.m_priv->addRef();
        m_naiveCost = ncmp.getNaiveCost();
        m_optimizedCost = ncmp.getOptimizedCost();
    }
}

bb::Matrix::Matrix()
    : m_priv(0)
    , m_ok(false)
{
    init(0, 0, "empty");
}

bb::Matrix::Matrix(int m, int n, std::string name)
    : m_priv(0)
    , m_ok(true)
{
    init(m, n, name);
}

void bb::Matrix::init(int m, int n, std::string name)
{
    m_priv = new MatrixPriv(m, n, name);
}

bb::Matrix::Matrix(const Matrix& m)
    : m_priv(0)
    , m_ok(false)
{
    if (m.m_priv) {
        m_priv = m.m_priv->addRef();
        m_ok = m.m_ok;
    }
}

bb::Matrix::~Matrix()
{
    if (m_priv)
        m_priv->release();
}

bool bb::Matrix::ok() const
{
    return m_ok;
}

std::string bb::Matrix::name() const
{
    return (m_priv ? m_priv->name : "BAD");
}

void bb::Matrix::setName(std::string name)
{
    if (m_priv)
        m_priv->name = name;
}

int bb::Matrix::rows() const
{
    return (m_priv ? m_priv->rows : 0);
}

int bb::Matrix::cols() const
{
    return (m_priv ? m_priv->cols : 0);
}

float* bb::Matrix::operator[](int r)
{
    assert(0 <= r && r < rows());
    return (m_priv ? m_priv->mat[r] : 0);
}

bb::Matrix& bb::Matrix::operator=(const Matrix& o)
{
    if (m_priv)
    {
        m_priv->release();
        m_priv = 0;
    }
    if (o.m_priv)
        m_priv = o.m_priv->addRef();
    return *this;
}

bb::Matrix bb::Matrix::mult(Matrix& b, size_t& numMults)
{
    Matrix result(rows(), b.cols(),
                  std::string("(") + name() + "*" + b.name() + ")");

    Matrix &a = *this;        // just for clean expressions below

    size_t mults = 0;

    for (int i=0; i < a.rows(); i++)
        for (int j=0; j < b.cols(); j++)
        {
            float sum = 0.0f;

            for (int k=0; k<a.cols(); k++)
                sum += a[i][k] * b[k][j];

            result[i][j] = sum;
            mults += a.cols();
        }

    numMults += mults;
    return result;
}

// ---------------- MatrixProduct

bool bb::MatrixProduct::calculate(Matrix &OUTresult)
{
    if (m_matrices.size() == 0)
        return false;

    if (m_matrices.size() == 1)
    {
        OUTresult = m_matrices[0];          // just copy it out
        return true;
    }

    const size_t last = m_matrices.size() - 1;

    // Confirm adjacent matrices are of same dimensionality.
    // OR, instead of bool return, throw an exception.
    for (size_t i=0; i<last; i++)
    {
        if (m_matrices[i].cols() != m_matrices[i+1].rows())
            return false;
    }

    m_numMults = 0;
    makeChainOrder(OUTresult);
    
    return true;
}

/*
 * Calculates the "cost" of matrix multiplication - the number of
 * multiplications two matrices require.
 */
size_t bb::MatrixProduct::getNaiveCost()
{
    size_t cost = 0;
    const size_t last = m_matrices.size() - 1;
    for (size_t i=0; i<last; i++)
    {
        // For two matrices mxn and nxp, the #mults is m*n*p
        cost += m_matrices[i].rows() * m_matrices[i].cols() * m_matrices[i+1].cols();
    }
    return cost;
}

size_t bb::MatrixProduct::getOptimizedCost()
{
    return m_numMults;
}

bb::MatrixProduct::MatrixProduct(const Matrix& m1, const Matrix& m2)
    : m_numMults(0)
{
    m_matrices.reserve(10);
    m_matrices.push_back(m1);
    m_matrices.push_back(m2);
}

bb::MatrixProduct& bb::MatrixProduct::operator*(const Matrix& m)
{
    m_matrices.push_back(m);
    return *this;
}

bb::MatrixProduct operator*(const bb::Matrix& m1, const bb::Matrix& m2)
{
    return bb::MatrixProduct(m1, m2);
}

/*
 * Utility functions to allocate and free 2D matrices of
 * integers.  Allocated arrays are initialized to zero.
 */
static int** alloc2DInt(uint r, uint c)
{
    int **result = new int*[r];
    for (uint i=0; i<r; i++)
    {
        result[i] = new int[c];
        ::memset(result[i], 0, sizeof(int) * c);
    }
    return result;
}

static void free2DInt(int **array, uint r)
{
    for (uint i=0; i<r; i++)
        delete[] array[i];
    delete[] array;
}

static void printMatrix(std::string msg, int **m, int rows, int cols)
{
    std::cout << "---- " << msg << std::endl;

    std::cout << "+-";
    for (int i=0; i<cols; i++)
        std::cout << "      ";
    std::cout << "-+\n";

    for (int r=0; r<rows; r++) {
        std::cout << "| ";
        for (int c=0; c<cols; c++)
            std::cout << std::setw(6) << m[r][c];
        std::cout << " |\n";
    }

    std::cout << "+-";
    for (int i=0; i<cols; i++)
        std::cout << "      ";
    std::cout << "-+\n";
}

/*
 * The following two functions come from
 * "Introduction to Algorithms", by Cormen, Leiserson, and Rivest, 1994.
 * p.306-308
 */
void bb::MatrixProduct::makeChainOrder(Matrix &OUTresult)
{
    size_t n = m_matrices.size();
    int *dims = new int[n + 1];

    // This function is called after the matrix chain has been
    // checked for correct dimensionality between adjacent matrices.
    // No need to check that [i-1].n == [i].m here.
    dims[0] = m_matrices[0].rows();
    for (int i=1; i<=n; i++)
        dims[i] = m_matrices[i-1].cols();

    // m[i,j] = Minimum number of scalar multiplications (i.e., cost)
    // needed to compute the matrix A[i]A[i+1]...A[j] = A[i..j]
    int **m = alloc2DInt(n+1, n+1);        // initialized to zero
    int **s = alloc2DInt(n+1, n+1);        // initialized to zero

    for (int len = 2; len <= n; len++)     // Subsequence lengths
    {
        for (int i = 1; i <= n - len + 1; i++)
        {
            int j = i + len - 1;
            m[i][j] = INT_MAX;
            for (int k = i; k <= j - 1; k++)
            {
                int cost = m[i][k] + m[k+1][j] + dims[i-1]*dims[k]*dims[j];
                if (cost < m[i][j])
                {
                    m[i][j] = cost;
                    s[i][j] = k; // Index of the subsequence split that achieved minimal cost
                }
            }
        }
    }

    // printMatrix("Cost Matrix", m, n+1, n+1);
    // printMatrix("Split Matrix", s, n+1, n+1);

    OUTresult = matrixChainMultiply(s, 1, n);

    free2DInt(s, n+1);
    free2DInt(m, n+1);
    delete[] dims;
}

bb::Matrix bb::MatrixProduct::matrixChainMultiply(int **s, int i, int j)
{
    if (j > i)
    {
        Matrix x = matrixChainMultiply(s, i, s[i][j]);
        Matrix y = matrixChainMultiply(s, s[i][j]+1, j);
        Matrix result = x.mult(y, m_numMults);
        return result;
    }
    else
        return m_matrices[i-1]; // [i] in book, [i-1] for 0-basis.
}

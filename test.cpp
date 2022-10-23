/*
 * Operand-collecting Matrix multiplication.
 * Test program.
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
#include "matrix.h"

void exampleChain()
{
    bb::Matrix a(30, 35, "A");
    bb::Matrix b(35, 15, "B");
    bb::Matrix c(15, 5,  "C");
    bb::Matrix d(5,  10, "D");
    bb::Matrix e(10, 20, "E");
    bb::Matrix f(20, 25, "F");

    std::cout << "exampleChain() start." << std::endl;

    bb::Matrix product = a * b * c * d * e * f;

    std::cout << "Final result is "
              << product.rows() << " x " << product.cols() << std::endl;
    std::cout << "Result name: " << product.name() << std::endl << std::endl;
    std::cout << "Simple Mult cost: " << product.getNaiveCost() << std::endl;
    std::cout << "Optimized Mult cost: " << product.getOptimizedCost() << std::endl;
    std::cout << "----------------\n";
}

void matrix4x4Comparison()
{
    bb::Matrix pt_in(4, 1, "A");  // col vector, a 3D point
    bb::Matrix to(4,4, "TO");     // translate to origin
    bb::Matrix rz(4,4, "RZ");     // rotate about Z
    bb::Matrix tb(4,4, "TB");     // translate back to world space

    std::cout << "matrix4x4Comparison() start." << std::endl;

    bb::Matrix pt_out = tb * rz * to * pt_in; // order of ops is right to left.

    std::cout << "Final result is "
              << pt_out.rows() << " x " << pt_out.cols() << std::endl;
    std::cout << "Result name: " << pt_out.name() << std::endl << std::endl;
    std::cout << "Simple Mult cost: " << pt_out.getNaiveCost() << std::endl;
    std::cout << "Optimized Mult cost: " << pt_out.getOptimizedCost() << std::endl;
}

int main(int argc, char **argv)
{
    exampleChain();
    matrix4x4Comparison();

    return 0;
}

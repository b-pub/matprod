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
    bb::Matrix a(30, 35, "A1");
    bb::Matrix b(35, 15, "A2");
    bb::Matrix c(15, 5,  "A3");
    bb::Matrix d(5,  10, "A4");
    bb::Matrix e(10, 20, "A5");
    bb::Matrix f(20, 25, "A6");

    std::cout << "exampleChain() start." << std::endl;

    bb::Matrix product = a * b * c * d * e * f;

    std::cout << "Final result is "
              << product.rows() << " x " << product.cols() << std::endl;
    std::cout << "Result name: " << product.name() << std::endl << std::endl;
}

int main(int argc, char **argv)
{
    exampleChain();

    return 0;
}

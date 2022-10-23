# matprod

An experiment in operand-collecting for matrix chain multiplication

The purpose of this code was to explore how a C++ expression
could be altered with operator overloading to collect operands
for later evaluation. This would allow appropriate optimizations
to be determined once all the operands are known. One type of
calculation that this is directly applicable to is matrix chain
multiplication, the application in this project.

# License

This software is released under the MIT license. See the LICENSE file
for details.

# Building

This project uses CMake 2.8 or later to generate the build files.
No other libraries are needed.

Once the repo is cloned, cd into the toplevel directory. When
using cmake, I like to locate the build directory in the toplevel
source directory. Create the build with:
```
    cd matprod
    mkdir build ; cd build
    cmake ..
```

CMake will then create build project files appropriate to the
platform.  By default on OS X or Linux, this creates
Makefiles. On Windows, it will create MSVC solution files. See
CMake's -G option to specify the desired type you need.

Once the build files are created, build it. On UNIX:
```
    make
```

The only executable defined in this project is a simple
test program called "mattest", which prints out the optimized
chain multiplication as a parenthesized expression.

Code written Oct 2016.
Brent Burton

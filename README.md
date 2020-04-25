## Finite Element Solver written in C++

This is a simple solver for the finite element method for simulating linear elastic deformation in 2D. It is written in C++ using [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page).

The core code is largely based on [Writing a FEM solver in less than 180 lines of code](https://podgorskiy.com/spblog/304/writing-a-fem-solver-in-less-the-180-lines-of-code), as great tutorial by Stanislav Pidhorskyi.

However, I expanded on it by adjusting the method for use in real-time applications and adding single-step methods for integrating the equations of motion.

### Building

The library can be added to your project as a submodule or by placing the repos contents in your project tree. It includes a `CMakeLists.txt` that defines the target `fem_solver` so you can include it in your CMake build by using `add_subdirectory`.

### Usage

TODO: Add some documentation.

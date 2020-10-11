# oml (Object Matrix Library)
This is a somewhat antiquated C++ matrix class which supports all the usual stuff:
- Matrix and Vector classes dynamic memory allocation.
- Fortran subscripting so we can easily pass matrices and vectors to FORTRAN code (lapack, arpack etc.).
- Shallow copy with copy on write symantics so passing matrices by value is *cheap* but safe.
- Expression templates for all sorts of overloated operators. As a result most of the code is in headers.
- Ascii and binary IO.
- Some linear algebra, SVD, Solve A*x=B, Eigen systems.
- Two and four tap XOR random number generators.
- Heap and shell sort.

This was used in the 1990s for electronic structure, monte carlo and electrochemical simulations coded in C++. As such it uses C++ 98. The expression templates could probably be greatly simplified with modern C++, but there are most likely much better matrix classes available now.  I still use it simply becuase I know how it works and it runs near machine speed.  I am currently using it to learn about Matrix Product States for solving 1D quantum spin systems.

The unit test files show how everything is used.
There some DOxygen comments throughout, but I havent tried to build the docs recently.
Right now it builds using the Code::Blocks IDE

# oml
This is a somewhat antiquated C++ matrix class which supports all the usuall stuff:
-Matrix and Vector classes dynamic memory allocation
-Fortran subscripting so we can easily pass matrices and vectors to FORTRAN code (lapack, arpack etc.)
-Shallow copy with copy on write symantics so passing matrices by values is *cheap* but safe
-Expression templates for all sorts of overloated operators
-Ascii and binary IO
-Somelinear algebra, SVD, Solve A*x=B, Eigen systems
-Two and four tap XOR random number generators
-Heap and sh ell sort

This was used 1990s for electronic structure, monte carlo and electrochemical simulations coded in C++. As such it uses C++ 98. The expression templates could probably be greatly simplified with modern C++, but there are most likely much better matrix classes available now.  I still use simply becuase I know how it works and it runs at machine speed.

The unit test files show how everything is used.

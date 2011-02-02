# oml/Makefile.am

PACKAGE = oml
VERSION = 1.3
pkgincludedir = $(includedir)/$(PACKAGE)
OPT=-O2

AM_CPPFLAGS=-I$(top_srcdir) -I$(top_builddir)
AUTOMAKE_OPTIONS=nostdinc

EXTRA_DIST=matsub.ci

#
# Names of all C++ source files, split in suitable sub categories.
#
SOURCES= src/matlimit.cc src/ran250.cc src/stream.cc src/veclimit.cc
DSOURCES= src/matlimit.cc src/ran250.cc src/stream.cc src/veclimit.cc
T_SOURCES= src/dmatrix.cc src/matrix.cc src/smatrix.cc 
N_SOURCES= numeric/cholsky.cc numeric/eigsort.cc numeric/gaussj.cc numeric/invmat.cc \
numeric/invtri.cc numeric/lubksb.cc numeric/lubksbm.cc numeric/ludcmp.cc numeric/lusolver.cc \
numeric/qldecomp.cc numeric/solvetri.cc numeric/svbksb.cc numeric/svdcmp.cc numeric/svsolver.cc \
numeric/tridiag.cc
CN_SOURCES= cnumeric/ch.cc cnumeric/epsilon.cc cnumeric/htribk.cc cnumeric/htridi.cc cnumeric/tql2.cc cnumeric/tqlrat.cc

#
#  Calculate names of object files that libtool creates.
#
T_OBJECTS= \
$(patsubst src/%.cc,libdouble_la-%.lo,$(T_SOURCES)) \
$(patsubst src/%.cc,libcomplex_la-%.lo,$(T_SOURCES)) \
$(patsubst numeric/%.cc,libdouble_la-%.lo,$(N_SOURCES)) \
$(patsubst cnumeric/%.cc,libcomplex_la-%.lo,$(CN_SOURCES))

TD_OBJECTS= \
$(patsubst src/%.cc,libdouble_d_la-%.lo,$(T_SOURCES)) \
$(patsubst src/%.cc,libcomplex_d_la-%.lo,$(T_SOURCES)) \
$(patsubst numeric/%.cc,libdouble_d_la-%.lo,$(N_SOURCES)) \
$(patsubst cnumeric/%.cc,libcomplex_d_la-%.lo,$(CN_SOURCES))

# 
# Names of library targets to build
#
lib_LTLIBRARIES = ../liboml.la ../libomld.la
#
# These are dummy libraries, just to fool libtool into making distinct object files.
#
EXTRA_LTLIBRARIES = libdouble.la libcomplex.la libdouble_d.la libcomplex_d.la

#
#  Non debug high performance versions.
#
___liboml_la_SOURCES = $(SOURCES)
___liboml_la_DEPENDENCIES=libdouble.la libcomplex.la $(T_OBJECTS)
___liboml_la_CXXFLAGS=$(OPT) -ftemplate-depth-30 -DNDEBUG
___liboml_la_LDFLAGS=-release $(VERSION)
___liboml_la_LIBADD  = $(T_OBJECTS)

libdouble_la_SOURCES = $(T_SOURCES) $(N_SOURCES)
libdouble_la_CXXFLAGS= -DTYPE=double $(OPT) -DNDEBUG

libcomplex_la_SOURCES = $(T_SOURCES) $(CN_SOURCES)
libcomplex_la_CXXFLAGS= -DTYPE='std::complex<double>' -include $(srcdir)/omlcomplex.h $(OPT) -DNDEBUG

#
#  Debug versions, full of asserts.
#
___libomld_la_SOURCES = $(DSOURCES)
___libomld_la_DEPENDENCIES=libdouble_d.la libcomplex_d.la $(TD_OBJECTS)
___libomld_la_CXXFLAGS=$(OPT) -ftemplate-depth-30 -DDEBUG
___libomld_la_LDFLAGS=-release $(VERSION)
___libomld_la_LIBADD  = $(TD_OBJECTS)

libdouble_d_la_SOURCES = $(T_SOURCES) $(N_SOURCES)
libdouble_d_la_CXXFLAGS= -DTYPE=double $(OPT) -DDEBUG

libcomplex_d_la_SOURCES = $(T_SOURCES) $(CN_SOURCES)
libcomplex_d_la_CXXFLAGS= -DTYPE='std::complex<double>' -include $(srcdir)/omlcomplex.h  $(OPT) -DDEBUG

nobase_pkginclude_HEADERS = array.h iterable.h mixtypes.h stopw.h array_io.h iterable_io.h numeric.h stream.h arrindex.h list.h  omlcomplex.h svsolver.h binio.h list_io.h tstream.h binop.h lusolver.h unop.h cnumeric.h matindex.h ran250.h unpickle.h cow.h matlimit.h random.h vecindex.h dmatrix.h matrix3d.h veclimit.h heapsort.h matrixalg.h rowcol.h vector3d.h histogram.h matrixbase.h shape.h vector.h indexable.h matrix.h shellsort.h vector_io.h indexsort.h matsub.h indext.h memops.h smatindex.h xpr.h io3d.h minmax.h smatrix.h subarray.h subvector.h submatrix.h

check_PROGRAMS = UT-Array-int UT-Array-double UT-Array-complex \
                 UT-List-int UT-List-double UT-List-complex \
                 UT-Vector-int UT-Vector-double UT-Vector-complex \
                 UT-Matrix-double UT-Matrix-complex UT-DMatrix-double UT-SMatrix-double \
                 UT-Vector3D-double UT-Matrix3D-double \
                 UT-Eigen-double UT-Random-double UT-SolveTriDiagonal


TESTS = $(check_PROGRAMS)
TESTS_CXXFLAGS=$(OPT)  
TESTS_LDADD   =-L$(top_builddir) -lomld
TESTS_CLEAN   =temp.dat tempb.dat ran.dat raw.dat

UT_Array_int_SOURCES       = UnitTests/UT-Array-int.cc UnitTests/UnitTest.cc UnitTests/UnitTest.H
UT_Array_int_LDADD         = $(TESTS_LDADD)
UT_Array_int_CXXFLAGS      = $(TESTS_CXXFLAGS)
UT_Array_double_SOURCES    = UnitTests/UT-Array-double.cc UnitTests/UnitTest.cc
UT_Array_double_LDADD      = $(TESTS_LDADD)
UT_Array_double_CXXFLAGS   = $(TESTS_CXXFLAGS)
UT_Array_complex_SOURCES   = UnitTests/UT-Array-complex.cc UnitTests/UnitTest.cc
UT_Array_complex_LDADD     = $(TESTS_LDADD)
UT_Array_complex_CXXFLAGS  = $(TESTS_CXXFLAGS)

UT_List_int_SOURCES        = UnitTests/UT-List-int.cc UnitTests/UnitTest.cc
UT_List_int_LDADD          = $(TESTS_LDADD)
UT_List_int_CXXFLAGS       = $(TESTS_CXXFLAGS)
UT_List_double_SOURCES     = UnitTests/UT-List-double.cc UnitTests/UnitTest.cc
UT_List_double_LDADD       = $(TESTS_LDADD)
UT_List_double_CXXFLAGS    = $(TESTS_CXXFLAGS)
UT_List_complex_SOURCES    = UnitTests/UT-List-complex.cc UnitTests/UnitTest.cc
UT_List_complex_LDADD      = $(TESTS_LDADD)
UT_List_complex_CXXFLAGS   = $(TESTS_CXXFLAGS)

UT_Vector_int_SOURCES      = UnitTests/UT-Vector-int.cc UnitTests/UnitTest.cc
UT_Vector_int_LDADD        = $(TESTS_LDADD)
UT_Vector_int_CXXFLAGS     = $(TESTS_CXXFLAGS)
UT_Vector_double_SOURCES   = UnitTests/UT-Vector-double.cc UnitTests/UnitTest.cc
UT_Vector_double_LDADD     = $(TESTS_LDADD)
UT_Vector_double_CXXFLAGS  = $(TESTS_CXXFLAGS)
UT_Vector_complex_SOURCES  = UnitTests/UT-Vector-complex.cc UnitTests/UnitTest.cc
UT_Vector_complex_LDADD    = $(TESTS_LDADD)
UT_Vector_complex_CXXFLAGS = $(TESTS_CXXFLAGS)

UT_Matrix_double_SOURCES   = UnitTests/UT-Matrix-double.cc UnitTests/UnitTest.cc
UT_Matrix_double_LDADD     = $(TESTS_LDADD)
UT_Matrix_double_CXXFLAGS  = $(TESTS_CXXFLAGS)
UT_DMatrix_double_SOURCES  = UnitTests/UT-DMatrix-double.cc UnitTests/UnitTest.cc
UT_DMatrix_double_LDADD    = $(TESTS_LDADD)
UT_DMatrix_double_CXXFLAGS = $(TESTS_CXXFLAGS)
UT_SMatrix_double_SOURCES  = UnitTests/UT-SMatrix-double.cc UnitTests/UnitTest.cc
UT_SMatrix_double_LDADD    = $(TESTS_LDADD)
UT_SMatrix_double_CXXFLAGS = $(TESTS_CXXFLAGS)
UT_Matrix_complex_SOURCES  = UnitTests/UT-Matrix-complex.cc UnitTests/UnitTest.cc
UT_Matrix_complex_LDADD    = $(TESTS_LDADD)
UT_Matrix_complex_CXXFLAGS = $(TESTS_CXXFLAGS)

UT_Vector3D_double_SOURCES = UnitTests/UT-Vector3D-double.cc UnitTests/UnitTest.cc
UT_Vector3D_double_LDADD   = $(TESTS_LDADD)
UT_Vector3D_double_CXXFLAGS= $(TESTS_CXXFLAGS)
UT_Matrix3D_double_SOURCES = UnitTests/UT-Matrix3D-double.cc UnitTests/UnitTest.cc
UT_Matrix3D_double_LDADD   = $(TESTS_LDADD)
UT_Matrix3D_double_CXXFLAGS= $(TESTS_CXXFLAGS)

UT_Eigen_double_SOURCES     = UnitTests/UT-Eigen-double.cc UnitTests/UnitTest.cc
UT_Eigen_double_LDADD       = $(TESTS_LDADD)
UT_Eigen_double_CXXFLAGS    = $(TESTS_CXXFLAGS)
UT_Random_double_SOURCES    = UnitTests/UT-Random-double.cc UnitTests/UnitTest.cc
UT_Random_double_LDADD      = $(TESTS_LDADD)
UT_Random_double_CXXFLAGS   = $(TESTS_CXXFLAGS)
UT_SolveTriDiagonal_SOURCES = UnitTests/UT-SolveTriDiagonal.cc UnitTests/UnitTest.cc
UT_SolveTriDiagonal_LDADD   = $(TESTS_LDADD)
UT_SolveTriDiagonal_CXXFLAGS= $(TESTS_CXXFLAGS)




CLEANFILES    =*~ UnitTests/*~ src/*~ numeric/*~ cnumeric/*~ $(EXTRA_LTLIBRARIES) $(TESTS_CLEAN)
DISTCLEANFILES=   $(EXTRA_LTLIBRARIES) $(TESTS_CLEAN) 
DOC_DIRS = .
DOC_EXAMPLES =
include $(top_srcdir)/Makefile.cmn





#include "oml/numeric/LapackCholsky.H"
#include "oml/smatrix.h"
#include "oml/matrix.h"


//
// https://netlib.org/lapack//explore-html/d4/d81/group__pptrf_ga767274ca963e0c56fcbbcd9ee23c4fc9.html#ga767274ca963e0c56fcbbcd9ee23c4fc9
// for detailed docs
// you also need to add -llapack to the link command

extern"C" {
void dpptrf_(char* uplo,int* N,double* A,int* INFO);
}

namespace oml {

Matrix<double> LapackCholsky(const SMatrix<double>& A)
{
    
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(!isnan(A));
    int N=A.GetNumRows(),info=0;
    char ul='L';
    SMatrix<double> V(A);
    dpptrf_(&ul,&N,&V(1,1),&info);
    if (info!=0)
        std::cerr << "Warning: LapackCholsky info=" << info << std::endl;
    assert(info==0);
    Matrix<double> Vfull(V);
    for (int i=1;i<=N;i++)
        for (int j=1;j<i;j++)
            Vfull(i,j)=0.0;
    return Vfull;
}


};

#include "LapackQRSolver.H"
#include "oml/matrix.h"
#include "oml/vector.h"
#include <complex>
//
// See http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html
// for detailed docs
// you also need to add -llapack to the link command
typedef std::complex<double> dcmplx;

extern"C" {
void dgeqrf_(int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void dgerqf_(int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void zgeqrf_(int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);
void zgerqf_(int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);

void dgeqlf_(int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void dgelqf_(int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void zgeqlf_(int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);
void zgelqf_(int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);

void dorgqr_(int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void dorgrq_(int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void zungqr_(int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);
void zungrq_(int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);

void dorgql_(int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void dorglq_(int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO);
void zungql_(int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);
void zunglq_(int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO);

}

template <class T> void xgeqrf  (int* M,int* N,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xgeqrf<double> (int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dgeqrf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xgeqrf<dcmplx> (int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zgeqrf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}

template <class T> void xgeqlf  (int* M,int* N,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xgeqlf<double> (int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dgeqlf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xgeqlf<dcmplx> (int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zgeqlf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}

template <class T> void xgelqf  (int* M,int* N,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xgelqf<double> (int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dgelqf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xgelqf<dcmplx> (int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zgelqf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}


template <class T> void xgerqf  (int* M,int* N,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xgerqf<double> (int* M,int* N,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dgerqf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xgerqf<dcmplx> (int* M,int* N,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zgerqf_(M,N,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}

template <class T> void xungqr  (int* M,int* N,int*	K,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xungqr<double> (int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dorgqr_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xungqr<dcmplx> (int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zungqr_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}

template <class T> void xungql  (int* M,int* N,int*	K,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xungql<double> (int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dorgql_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xungql<dcmplx> (int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zungql_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}

template <class T> void xungrq  (int* M,int* N,int*	K,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xungrq<double> (int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dorgrq_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xungrq<dcmplx> (int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zungrq_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}

template <class T> void xunglq  (int* M,int* N,int*	K,T* Q,int* LDA,T* TAU,T* WORK,int* LWORK,int* INFO);
template <> void xunglq<double> (int* M,int* N,int*	K,double* Q,int* LDA,double* TAU,double* WORK,int* LWORK,int* INFO)
{
    dorglq_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //double
}
template <> void xunglq<dcmplx> (int* M,int* N,int*	K,dcmplx* Q,int* LDA,dcmplx* TAU,dcmplx* WORK,int* LWORK,int* INFO)
{
    zunglq_(M,N,K,Q,LDA,TAU,WORK,LWORK,INFO); //complex<double>
}


template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::SolveThinQR(const MatrixT& Ain)
{
    assert(Ain.GetLimits().Row.Low==1);
    assert(Ain.GetLimits().Col.Low==1);
    int M=Ain.GetNumRows(),N=Ain.GetNumCols(),mn=Min(M,N);
    int info=0,lwork=-1;
    Vector<T> tau(mn),work(1);
    Matrix<T> Q(Ain);
    //
    //  Initial call to see how much work space is needed
    //
    xgeqrf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    //
    //  Now do the actual QR work
    //
    xgeqrf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    //
    //  Grab R before xungqr clobbers it.
    //
    Matrix<T> R(mn,N);
    for (int i=1;i<=mn;i++)
    {
        for (int j=1;j<i;j++)
            R(i,j)=0.0;
        for (int j=i;j<=N;j++)
            R(i,j)=Q(i,j);

    }
    //
    //  unpack Q
    //
    lwork=-1;
    info=0;
    xungqr<T>(&M,&mn,&mn,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    xungqr<T>(&M,&mn,&mn,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    if (mn<N) Q.SetLimits(M,mn,true);
    return std::make_tuple(std::move(Q),std::move(R)); //return [Q,R]
}

template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::SolveThinQL(const MatrixT& Ain)
{
    assert(Ain.GetLimits().Row.Low==1);
    assert(Ain.GetLimits().Col.Low==1);
    int M=Ain.GetNumRows(),N=Ain.GetNumCols(),mn=Min(M,N);
    int info=0,lwork=-1;
    Vector<T> tau(mn),work(1);
    Matrix<T> Q(Ain);
    //
    //  Initial call to see how much work space is needed
    //
    xgeqlf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    //
    //  Now do the actual QL work
    //
    xgeqlf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    //
    //  Grab L before xungql clobbers it.
    //
    Matrix<T> L(mn,N);
    for (int i=1;i<=mn;i++)
    {
        for (int j=1;j<=i+N-mn;j++)
            L(i,j)=Q(i+M-mn,j);
        for (int j=i+1+N-mn;j<=N;j++)
            L(i,j)=0.0;
    }
    //
    // Now we need shift the orth vectors from the right side of Q over the left side, before
    // shrinking Q and extracting the reflectors.
    //
    if (mn<N)
    { //
      //
        for (int i=1;i<=M;i++)
            for (int j=1;j<=mn;j++)
                Q(i,j)=Q(i,j+N-mn);
        Q.SetLimits(M,mn,true);
    }
    //
    //  unpack Q  if M<N then the reflectors on the right side of Q, so we need
    //  a shifted Q to xungql
    //
    lwork=-1;
    info=0;
    xungql<T>(&M,&mn,&mn,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    xungql<T>(&M,&mn,&mn,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);

    return std::make_tuple(std::move(Q),std::move(L)); //return [Q,R]
}

template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::SolveThinRQ(const MatrixT& Ain)
{
    assert(Ain.GetLimits().Row.Low==1);
    assert(Ain.GetLimits().Col.Low==1);
    int M=Ain.GetNumRows(),N=Ain.GetNumCols(),mn=Min(M,N);
    int info=0,lwork=-1;
    Vector<T> tau(mn),work(1);
    Matrix<T> Q(Ain);
    //
    //  Initial call to see how much work space is needed
    //
    xgerqf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    //
    //  Now do the actual QR work
    //
    xgerqf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    //
    //  Grab R before xungqr clobbers it.
    //
    Matrix<T> R(M,mn);
    for (int j=1;j<=mn;j++)
    {
        for (int i=1;i<=j+M-mn;i++)
            R(i,j)=Q(i,j+N-mn);
        for (int i=j+1+M-mn;i<=M;i++)
            R(i,j)=0.0;
    }
    //
    // If M>N we need shift the orth vectors from the bottom of Q up to top before
    // unpacking the reflectors.
    //
    if (mn<M)
    {
        for (int j=1;j<=N;j++)
            for (int i=1;i<=mn;i++)
                Q(i,j)=Q(i+M-mn,j);
        Q.SetLimits(mn,N,true);
    }
    //
    //  unpack Q
    //
    lwork=-1;
    info=0;
    xungrq<T>(&mn,&N,&mn,&Q(1,1),&mn,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    xungrq<T>(&mn,&N,&mn,&Q(1,1),&mn,&tau(1),&work(1),&lwork,&info);
    assert(info==0);

    return std::make_tuple(std::move(R),std::move(Q)); //Return [R,Q]
}

template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::SolveThinLQ(const MatrixT& Ain)
{
    assert(Ain.GetLimits().Row.Low==1);
    assert(Ain.GetLimits().Col.Low==1);
    int M=Ain.GetNumRows(),N=Ain.GetNumCols(),mn=Min(M,N);
    int info=0,lwork=-1;
    Vector<T> tau(mn),work(1);
    Matrix<T> Q(Ain);
    //
    //  Initial call to see how much work space is needed
    //
    xgelqf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    //
    //  Now do the actual QR work
    //
    xgelqf<T>(&M,&N,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    //
    //  Grab R before xungqr clobbers it.
    //
    Matrix<T> L(M,mn);
    for (int j=1;j<=mn;j++)
    {
        for (int i=1;i<j;i++)
            L(i,j)=0.0;
        for (int i=j;i<=M;i++)
            L(i,j)=Q(i,j);
    }
    //
    //  unpack Q
    //
    lwork=-1;
    info=0;
    xunglq<T>(&mn,&N,&mn,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    lwork=static_cast<int>(real(work(1)));
    work.SetLimits(lwork);
    xunglq<T>(&mn,&N,&mn,&Q(1,1),&M,&tau(1),&work(1),&lwork,&info);
    assert(info==0);
    if (mn<M) Q.SetLimits(mn,N,true);
    return std::make_tuple(std::move(L),std::move(Q)); //Return [L,Q]
}

//
//  Hunt for the zero pivots in L and remove those rows from L and columns from Q
//
template <class T> bool ShrinkQR(Matrix<T>& Q, Matrix<T>& R,double eps)
{
    Vector<T> diag=R.GetDiagonal();
    std::vector<index_t> remove;
    for (index_t i:diag.indices())
    {
        if (fabs(diag(i))<=eps) //fast feasibility study
        {
            double s=Sum(fabs(R.GetRow(i))); //Now check the whole row
            if (s<=eps) remove.push_back(i);
        }
    }
    for (auto i=remove.rbegin();i!=remove.rend();i++)
    {
        assert(Sum(fabs(R.GetRow(*i)))<eps); //Make sure we got the correct row
        R.RemoveRow   (*i);
        Q.RemoveColumn(*i);
    }
    return remove.size()>0;
}
//
//  Hunt for the zero pivots in L and remove those cols from L and rows from Q
//
template <class T> bool ShrinkRQ(Matrix<T>& R, Matrix<T>& Q,double eps)
{
    Vector<T> diag=R.GetDiagonal();
    std::vector<index_t> remove;
    for (index_t i:diag.indices())
    {
        if (fabs(diag(i))<=eps) //fast feasibility study
        {
            double s=Sum(fabs(R.GetColumn(i))); //Now check the whole col
            if (s<=eps) remove.push_back(i);
        }
    }
    for (auto i=remove.rbegin();i!=remove.rend();i++)
    {
        assert(Sum(fabs(R.GetColumn(*i)))<eps); //Make sure we got the correct row
        R.RemoveColumn(*i);
        Q.RemoveRow   (*i);
    }
    return remove.size()>0;
}

template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::
SolveRankRevealingQR(const MatrixT& A, double eps)
{
    auto [Q,R]=SolveThinQR(A);
    ShrinkQR(Q,R,eps);
    return std::make_tuple(Q,R);
}
template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::
SolveRankRevealingRQ(const MatrixT& A, double eps)
{
    auto [R,Q]=SolveThinRQ(A);
    ShrinkRQ(R,Q,eps);
    return std::make_tuple(R,Q);
}


template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::
SolveRankRevealingQL(const MatrixT& A, double eps)
{
    auto [Q,L]=SolveThinQL(A);
    ShrinkQR(Q,L,eps);
    return std::make_tuple(Q,L);
}

template <class T> typename LapackQRSolver<T>::QRType LapackQRSolver<T>::
SolveRankRevealingLQ(const MatrixT& A, double eps)
{
    auto [L,Q]=SolveThinLQ(A);
    ShrinkRQ(L,Q,eps);
    return std::make_tuple(L,Q);
}


template class LapackQRSolver<double>;
template class LapackQRSolver<dcmplx>;

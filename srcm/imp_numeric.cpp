module;
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>

module oml.NumericalRecipes:Imp;
import oml.Vector;
import oml.SMatrix;
import oml.Matrix;
import oml.IndexSort;

//-------------------------------------------------------------
//
//     estimate unit roundoff in quantities of size x.
//
//    this program should function properly on all systems
//    satisfying the following two assumptions,
//       1.  the base used in representing floating point
//           numbers is not a power of three.
//       2.  the quantity  a  in statement 10 is represented to 
//           the accuracy used in floating point variables
//           that are stored in memory.
//    the statement number 10 and the go to 10 are intended to
//    force optimizing compilers to generate code satisfying 
//    assumption 2.
//    under these assumptions, it should be true that,
//           a  is not exactly equal to four-thirds,
//           b  has a zero for its last bit or digit,
//           c  is not exactly equal to one,
//           eps  measures the separation of 1.0 from
//                the next larger floating point number.
//
double epsilon (double x)
{
  double a,b,c,eps;
  a = 4.0/3.0;
  do
    {
      b = a - 1.0;
      c = b + b + b;
      eps = fabs(c-1.0);
    } while (eps == 0.0);
	
  return eps*fabs(x);
}

//###########################################################################
//
//  Cholsky decomposition of a symmetric-positive definite matrix into
//  upper and lower triangular parts.  A -> U * ~U
//
template <class T> void Cholsky(Matrix<T>& A)
{
    assert(A.GetRowLimits()==A.GetColLimits());
    assert(IsSymmetric(A));
    assert(A.GetRowLow()==1);
    assert(A.GetColLow()==1);

    index_t n=A.GetNumRows();

    typename Matrix<T>::Subscriptor a(A);

    for(index_t j=n; j>=1; j--)
    {
        T temp=0.0;
        for(index_t k=j+1; k<=n; k++) temp+=a(j,k)*a(j,k);
        if (temp>a(j,j))
        {
            std::cerr << "Cholsky(SMatrix<T>& A): Matrix was not positive definite" << std::endl;
            assert(false);
        }
         a(j,j)=sqrt(a(j,j)-temp);
        for(index_t i=1; i<=j-1; i++)
        {
            temp=0.0;
            for(index_t k=j+1; k<=n; k++) temp+=a(i,k)*a(j,k);
            temp=a(i,j)-temp;
            if (temp!=0) a(i,j)=temp/a(j,j);
            else a(i,j)=0.0;
        }
    }
    for(index_t j=1; j<=n; j++)
        for(index_t i=j+1; i<=n; i++) a(i,j)=0.0;
}

//###########################################################################
//
//  Cholsky decomposition of a symmetric-positive definite matrix into
//  upper and lower triangular parts.  A -> U * ~U
//  Symmetric version, works on upper part, lower part will not be 0.0's.
//
//#include "oml/imp/isnan.h"
#define UPPER_ONLY
template <class T> void Cholsky(SMatrix<T>& A)
{
    assert(A.GetRowLow()==1);
    assert(A.GetColLow()==1);

    index_t n=A.GetNumRows();

    typename SMatrix<T>::Subscriptor a(A);

    for(index_t j=n; j>=1; j--)
    {
        T temp=0.0;
        for(index_t k=j+1; k<=n; k++) temp+=a(j,k)*a(j,k);
        if (temp>a(j,j))
        {
            std::cerr << "Cholsky(SMatrix<T>& A): Matrix was not positive definite" << std::endl;
            assert(false);
        }
        a(j,j)=sqrt(a(j,j)-temp);
        for(index_t i=1; i<=j-1; i++)
        {
            temp=0.0;
            for(index_t k=j+1; k<=n; k++) temp+=a(i,k)*a(j,k);
            temp=a(i,j)-temp;
            if (temp!=0) a(i,j)=temp/a(j,j);
            else a(i,j)=0.0;
        }
    }
}
#undef UPPER_ONLY
//##########################################################################
//
//  Sorts eigen values in ascending order, vectors are correspondingly
//  rearranged
//

template <class T, class M> void EigenSort(M& A, Vector<T>& EigenValues)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  typename M        ::Subscriptor a(A);
  typename Vector<T>::Subscriptor e(EigenValues);

  index_t n=A.GetNumRows(), k;

  for (index_t i=1;i<n;i++)
  {
    T p=e(k=i);
    for (index_t j=i+1;j<=n;j++) if (e(j) <= p) p=e(k=j);
    if (k != i)
    {
      e(k)=e(i); e(i)=p;
      for (index_t j=1;j<=n;j++) {p=a(j,i); a(j,i)=a(j,k); a(j,k)=p;}
    }
  }
}

//----------------------------------------------------------------------------
//
//  Matrix inversion.
//
template <class T> Matrix<T> InvertSquare(const Matrix<T>& m)
{
  Matrix<T> temp(m);
  std::vector<index_t> SwapIndex(temp.GetNumRows());
  Matrix<T> inv(temp.GetLimits());
  Unit(inv);
  T sign;
  LUDecomp (temp,SwapIndex,sign);
  LUBackSub(temp,inv,SwapIndex);
  return inv;
}


template <class T> Matrix<T> InvertSymmetric(const Matrix<T>& m)
{
  Matrix<T> ret(m);
  Vector<T> normal=T(1)/sqrt(m.GetDiagonal());
  //m=DirectMultiply(m,OuterProduct(normal));
//  Vector<T> normal=Normalize(ret); //Normalize.
  Cholsky(ret);             //Decompose.
  InvertTriangular(ret);    //Invert upper triangular matrix.
  ret=Matrix<T>(~ret*ret);             //Recompose.
  ret=DirectMultiply(ret,OuterProduct(normal)); //Unormalize.
  //Normalize(ret,normal);      
  return ret;
}

#define UPPER_ONLY
template <class T> SMatrix<T> InvertSymmetric(const SMatrix<T>& m)
{
  int n=m.GetNumRows();
  SMatrix<T> ret (m.GetLimits());
  Fill(ret,0.0);

  SMatrix<T> temp(m);
  //Vector<T> norm=Normalize(temp); //Normalize.
  Vector<T> norm=T(1)/sqrt(temp.GetDiagonal());
  temp=DirectMultiply(temp,OuterProduct(norm));
  
  Cholsky(temp);             //Decompose.
  InvertTriangular(temp);    //Invert upper triangular matrix.
  typename SMatrix<T>::Subscriptor      sr(ret);

  for (int i=1;i<=n;i++)
    for (int j=i;j<=n;j++)
      for (int k=1;k<=i;k++) sr(i,j)+=temp(k,i)*temp(k,j);

  //Normalize(ret,norm);      //Unormalize.
  ret=DirectMultiply(ret,OuterProduct(norm)); //Unormalize.
  
  return ret;
}
#undef _UPPER_ONLY
//###########################################################################
//
//  Invert an upper triangular matrix.
//
//
template <class T> void InvertTriangular(Matrix<T>& A)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();

  typename Matrix<T>::Subscriptor a(A);

  for(index_t i=1; i<=n; i++) a(i,i)=1.0/a(i,i);
  for(index_t j=n; j>=2; j--)
  {
    for (index_t i=j-1; i>=1; i--)
    {
      T temp=0.0;
      for(index_t k=i+1; k<=j; k++) temp+=a(i,k)*a(k,j);
      a(i,j)=-a(i,i)*temp;
    }
  }
}

//###########################################################################
//
//  Invert a symmetrically stored upper triangular matrix.
//
//
#define UPPER_ONLY
template <class T> void InvertTriangular(SMatrix<T>& A)
{
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();

  typename SMatrix<T>::Subscriptor a(A);

  for(index_t i=1; i<=n; i++) a(i,i)=1.0/a(i,i);
  for(index_t j=n; j>=2; j--)
  {
    for (index_t i=j-1; i>=1; i--)
    {
      T temp=0.0;
      for(index_t k=i+1; k<=j; k++) temp+=a(i,k)*a(k,j);
      a(i,j)=-a(i,i)*temp;
    }
  }
}
#undef UPPER_ONLY
template <class T> void LUBackSub(const Matrix<T>& a, Vector<T>& B,const std::vector<index_t>& Index)
{
  assert(a.GetRowLimits()==a.GetColLimits());
  assert(B.GetLimits   ()==a.GetRowLimits());
  assert(a.GetRowLow()==1);
  assert(a.GetColLow()==1);
  int i,ii=0,ip,j;
  T sum;

  typename Vector<T>::Subscriptor      b(B);
  index_t n=B.size();

  for (i=1;i<=n;i++)
  {
    ip=Index[i-1];
    sum=b(ip);
    b(ip)=b(i);
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a(i,j)*b(j);
    else
      if (sum) ii=i;
    b(i)=sum;
  }
  for (i=n;i>=1;i--)
  {
    sum=b(i);
    for (j=i+1;j<=n;j++) sum -= a(i,j)*b(j);
    b(i)=sum/a(i,i);
  }
}

//###########################################################################
//
//  LU backsubstitution of a whole matrix.
//

template <class T> void LUBackSub(const Matrix<T>& a, Matrix<T>& B, const std::vector<index_t>& Index)
{
  assert(a.GetRowLimits()==a.GetColLimits());
  assert(a.GetColLimits()==B.GetRowLimits());
  assert(a.GetRowLow()==1);
  assert(a.GetColLow()==1);
  assert(B.GetRowLow()==1);
  assert(B.GetColLow()==1);


  typename Matrix<T>::Subscriptor b(B);
  index_t n =a.GetNumRows();

  for (index_t isub=1; isub<=n; isub++)
  {
    int i,ii=0,ip,j;
    T sum;
    for (i=1; i<=n; i++)
    {
      ip=Index[i-1];
      sum=b(ip,isub);
      b(ip,isub)=b(i,isub);
      if (ii)
	for (j=ii; j<=i-1; j++) sum -= a(i,j)*b(j,isub);
      else
	if (sum) ii=i;
      b(i,isub)=sum;
    } //i
    for (i=n;i>=1;i--)
    {
      sum=b(i,isub);
      for (j=i+1; j<=n; j++) sum -= a(i,j)*b(j,isub);
      b(i,isub)=sum/a(i,i);
    } // i
  } // isub
}

template <class T> bool LUDecomp(Matrix<T>& A, std::vector<index_t>& ipiv ,T& d)
{
    assert(A.GetRowLimits()==A.GetColLimits());
    assert(A.GetRowLow()==1);
    assert(A.GetColLow()==1);
    T big,dum,sum,temp;
    index_t imax;
    const double TINY=1.0e-20;

    Vector<T> V(A.GetRowLimits());
    index_t n=V.size();

    typename Matrix<T>::Subscriptor a (A);
    typename Vector<T>::Subscriptor vv(V);

    d=1.0;
    for (index_t i=1; i<=n; i++)
    {
        big=0.0;
        for (index_t j=1; j<=n; j++)
            if ((temp=fabs(a(i,j))) > big) big=temp;
        if (big == 0.0)
        {
            std::cerr << "Singular matrix in routine LUDCMP" << std::endl;
            return false;
        }
        vv(i)=1.0/big;
    }
    for (index_t j=1; j<=n; j++)
    {
        for (index_t i=1; i<j; i++)
        {
            sum=a(i,j);
            for (index_t k=1; k<i; k++) sum -= a(i,k)*a(k,j);
            a(i,j)=sum;
        }
        big=-1.0;
        imax=0;
        for (index_t i=j; i<=n; i++)
        {
            sum=a(i,j);
            for (index_t k=1; k<j; k++)
                sum -= a(i,k)*a(k,j);
            a(i,j)=sum;
            if ( (dum=vv(i)*fabs(sum)) >= big)
            {
                big=dum;
                imax=i;
            }
        }
        assert(imax>0);
        if (j != imax)
        {
            for (index_t k=1; k<=n; k++)
            {
                dum=a(imax,k);
                a(imax,k)=a(j,k);
                a(j,k)=dum;
            }
            d = -d;
            vv(imax)=vv(j);
        }
        ipiv[j-1]=imax;
        if (a(j,j) == 0.0) a(j,j)=TINY;
        if (j != n)
        {
            dum=1.0/(a(j,j));
            for (index_t i=j+1; i<=n; i++) a(i,j) *= dum;
        }
    }
    return true;
}
//##########################################################################
//
//  QL algorithm
//  Calculates Eigenvalues and Eigenvectors of a tridiagonal matrix
//
//  Diagonal    - Input std::vector of diagonal elements
//              - Output std::vector of eigen values
//  OffDiagonal - Input std::vector of Off diagonal elements
//  M           - Output matrix containing eigen std::vectors as columns.
//
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))


template <class T, class M> void QLDecomp(M& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();

  T g,r,s,c,p,f,b;

  typename M        ::Subscriptor a (A);
  typename Vector<T>::Subscriptor d (   Diagonal);
  typename Vector<T>::Subscriptor od(OffDiagonal);

  for (index_t i=2;i<=n;i++) od(i-1)=od(i);
  od(n)=0.0;
  for (index_t l=1;l<=n;l++)
  {
    index_t iter=0;
    index_t m;
    do
    {
      for (m=l;m<=n-1;m++)
      {
	T dd=fabs(d(m))+fabs(d(m+1));
	if (fabs(od(m))+dd == dd) break;
      }
      if (m != l)
      {
	if (iter++ == 30)
        {
          std::cerr << "Too many iterations in QL_Decomp";
          assert(false);
        }
	g=(d(l+1)-d(l))/(2.0*od(l));
	r=sqrt((g*g)+1.0);
	g=d(m)-d(l)+od(l)/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (index_t i=m-1;i>=l;i--)
	{
	  f=s*od(i);
	  b=c*od(i);
	  if (fabs(f) >= fabs(g))
	  {
	    c=g/f;
	    r=sqrt((c*c)+1.0);
	    od(i+1)=f*r;
	    c *= (s=1.0/r);
	  }
	  else
	  {
	    s=f/g;
	    r=sqrt((s*s)+1.0);
	    od(i+1)=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d(i+1)-p;
	  r=(d(i)-g)*s+2.0*c*b;
	  p=s*r;
	  d(i+1)=g+p;
	  g=c*r-b;
	  for (index_t k=1;k<=n;k++)
	  {
	    f=a(k,i+1);
	    a(k,i+1)=s*a(k,i)+c*f;
	    a(k,i  )=c*a(k,i)-s*f;
	  }
	}
	d(l)=d(l)-p;
	od(l)=g;
	od(m)=0.0;
      }
    } while (m != l);
  }
}
#undef SIGN

template <class T> void SVBackSub(const Matrix<T>& u, const Vector<T>& w,const Matrix<T>& v, const Vector<T>& b, Vector<T>& X)
{
   assert(u.GetColLimits()==w.GetLimits());
   assert(u.GetColLimits()==X.GetLimits());
   assert(u.GetColLimits()==v.GetRowLimits());
   assert(u.GetColLimits()==v.GetColLimits());
   assert(u.GetRowLimits()==b.GetLimits());

   int jj,j,i;
   T    s;
   index_t m=u.GetNumRows();
   index_t n=u.GetNumCols();

   typename Vector<T>::Subscriptor x(X);

   Vector<T> tmp(n);

   for (j=1;j<=n;j++)
   {
     s=0.0;
     if (w(j))
     {
       for (i=1;i<=m;i++) s += u(i,j)*b(i);
       s /= w(j);
     }
     tmp(j)=s;
   }
   for (j=1;j<=n;j++)
   {
     s=0.0;
     for (jj=1;jj<=n;jj++) s += v(j,jj)*tmp(jj);
     x(j)=s;
   }
}

template <class T> T PYTHAG(T a,T b) 
{
  T at=fabs(a);
  T bt=fabs(b);
  T ct=at>bt ? bt/at : at/bt;
  return at>bt ? at*sqrt(1.0+ct*ct) : (bt ? bt*sqrt(1.0+ct*ct) : 0.0);
}
// static double at,bt,ct;
// #define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? 
// (ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))


template <class T> T MAX(T a,T b) {return a > b ? a :b;}
// static double maxarg1,maxarg2;
// #define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?
//    (maxarg1) : (maxarg2))


template <class T> T SIGN(T a, T b) {return b >= 0.0 ? fabs(a) : -fabs(a);}
// #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

template <class T, class M> void SVDecomp(M& A, Vector<T>& W, M& V)
{
  assert(A.GetColLimits()==W.GetLimits());
  assert(A.GetColLimits()==V.GetRowLimits());
  assert(A.GetColLimits()==V.GetColLimits());

  int flag,l=0,nm=0;
  T c,f,h,s,x,y,z;
  T anorm=0.0,g=0.0,scale=0.0;
  index_t m=A.GetNumRows();
  index_t n=A.GetNumCols();

  if (m < n)
    {
      std::cerr << "SVDCMP: Resizing A with extra zero rows" << std::endl;
      A.SetLimits(n,n,true); //Add the requited rows automatically
      for (int i=m+1;i<=n;i++)
        for (int j=1;j<=n;j++)
            A(i,j)=T(0.0);
    }

  Vector<T> RV1(n);

  typename M::Subscriptor a(A);
  typename M::Subscriptor v(V);
  typename Vector<T>::Subscriptor w(W);
  typename Vector<T>::Subscriptor rv1(RV1);

  for (index_t i=1;i<=n;i++)
    {
      l=i+1;
      rv1(i)=scale*g;
      g=s=scale=0.0;
      if (i <= m)
	{
	  for (index_t k=i;k<=m;k++) scale += fabs(a(k,i));
	  if (scale)
	    {
	      for (index_t k=i;k<=m;k++)
		{
		  a(k,i) /= scale;
		  s += a(k,i)*a(k,i);
		}
	      f=a(i,i);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      a(i,i)=f-g;
	      if (i != n)
		{
		  for (index_t j=l;j<=n;j++)
		    {
		      s=0.0;
		      for (index_t k=i;k<=m;k++) s += a(k,i)*a(k,j);
		      f=s/h;
		      for (index_t k=i;k<=m;k++) a(k,j) += f*a(k,i);
		    }
		}
	      for (index_t k=i;k<=m;k++) a(k,i) *= scale;
	    }
	}
      w(i)=scale*g;
      g=s=scale=0.0;
      if (i <= m && i != n)
	{
	  for (index_t k=l;k<=n;k++) scale += fabs(a(i,k));
	  if (scale)
	    {
	      for (index_t k=l;k<=n;k++)
		{
		  a(i,k) /= scale;
		  s += a(i,k)*a(i,k);
		}
	      f=a(i,l);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      a(i,l)=f-g;
	      for (index_t k=l;k<=n;k++) rv1(k)=a(i,k)/h;
	      if (i != m) {
		for (index_t j=l;j<=m;j++)
		  {
		    s=0.0;
		    for (index_t k=l;k<=n;k++) s += a(j,k)*a(i,k);
		    for (index_t k=l;k<=n;k++) a(j,k) += s*rv1(k);
		  }
	      }
	      for (index_t k=l;k<=n;k++) a(i,k) *= scale;
	    }
	}
      anorm=MAX(anorm,(fabs(w(i))+fabs(rv1(i))));
    }
  for (index_t i=n;i>=1;i--)
    {
      if (i < n)
	{
	  if (g)
	    {
	      for (index_t j=l;j<=n;j++) v(j,i)=(a(i,j)/a(i,l))/g;
	      for (index_t j=l;j<=n;j++)
		{
		  s=0.0;
		  for (index_t k=l;k<=n;k++) s += a(i,k)*v(k,j);
		  for (index_t k=l;k<=n;k++) v(k,j) += s*v(k,i);
		}
	    }
	  for (index_t j=l;j<=n;j++) v(i,j)=v(j,i)=0.0;
	}
      v(i,i)=1.0;
      g=rv1(i);
      l=i;
    }
  for (index_t i=n;i>=1;i--)
    {
      l=i+1;
      g=w(i);
      if (i < n)
	for (index_t j=l;j<=n;j++) a(i,j)=0.0;
      if (g)
	{
	  g=1.0/g;
	  if (i != n) {
	    for (index_t j=l;j<=n;j++)
	      {
		s=0.0;
		for (index_t k=l;k<=m;k++) s += a(k,i)*a(k,j);
		f=(s/a(i,i))*g;
		for (index_t k=i;k<=m;k++) a(k,j) += f*a(k,i);
	      }
	  }
	  for (index_t j=i;j<=m;j++) a(j,i) *= g;
	}
      else
	{
	  for (index_t j=i;j<=m;j++) a(j,i)=0.0;
	}
      ++a(i,i);
    }
  for (index_t k=n;k>=1;k--)
    {
      for (index_t its=1;its<=30;its++)
	{
	  flag=1;
	  for (l=k;l>=1;l--)
	    {
	      nm=l-1;
	      if (fabs(rv1(l))+anorm == anorm)
		{
		  flag=0;
		  break;
		}
	      if (fabs(w(nm))+anorm == anorm) break;
	    }
	  if (flag)
	    {
	      c=0.0;
	      s=1.0;
	      for (index_t i=l;i<=k;i++)
		{
		  f=s*rv1(i);
		  if (fabs(f)+anorm != anorm)
		    {
		      g=w(i);
		      h=PYTHAG(f,g);
		      w(i)=h;
		      h=1.0/h;
		      c=g*h;
		      s=(-f*h);
		      for (index_t j=1;j<=m;j++)
			{
			  y=a(j,nm);
			  z=a(j,i);
			  a(j,nm)=y*c+z*s;
			  a(j,i)=z*c-y*s;
			}
		    }
		}
	    }
	  z=w(k);
	  if (l == k)
	    {
	      if (z < 0.0)
		{
		  w(k) = -z;
		  for (index_t j=1;j<=n;j++) v(j,k)=(-v(j,k));
		}
	      break;
	    }
	  if (its == 30)
	    {
	      std::cerr << "No convergence in 30 SVDCMP iterations" << std::endl;
	    }
	  x=w(l);
	  nm=k-1;
	  y=w(nm);
	  g=rv1(nm);
	  h=rv1(k);
	  f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g=PYTHAG(f,1.0);
	  f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	  c=s=1.0;
	  for (index_t j=l;j<=nm;j++)
	    {
	      index_t i=j+1;
	      g=rv1(i);
	      y=w(i);
	      h=s*g;
	      g=c*g;
	      z=PYTHAG(f,h);
	      rv1(j)=z;
	      c=f/z;
	      s=h/z;
	      f=x*c+g*s;
	      g=g*c-x*s;
	      h=y*s;
	      y=y*c;
	      for (index_t jj=1;jj<=n;jj++)
		{
		  x=v(jj,j);
		  z=v(jj,i);
		  v(jj,j)=x*c+z*s;
		  v(jj,i)=z*c-x*s;
		}
	      z=PYTHAG(f,h);
	      w(j)=z;
	      if (z)
		{
		  z=1.0/z;
		  c=f*z;
		  s=h*z;
		}
	      f=(c*g)+(s*y);
	      x=(c*y)-(s*g);
	      for (index_t jj=1;jj<=m;jj++)
		{
		  y=a(jj,j);
		  z=a(jj,i);
		  a(jj,j)=y*c+z*s;
		  a(jj,i)=z*c-y*s;
		}
	    }
	  rv1(l)=0.0;
	  rv1(k)=f;
	  w(k)=x;
	}
    }
    //
    // Sort for descending singular values
    //
    std::vector<index_t> index=MakeDescendingIndex(W); //order of s is preserved.
    W.ReIndex(index);
    A.ReIndexColumns(index);
    V.ReIndexColumns(index);
    //
    //  Reshape to remove temp junk
    //
    if (m < n)
    {
        W.SetLimits(static_cast<size_t>(m),true);
        A.SetLimits(m,m,true);
        V.SetLimits(m,n,true);
    }
}

#undef SIGN
#undef MAX
#undef PYTHAG


//###########################################################################
//
//  Housholder reduction of a real symmetric matrix, to tri-diagonal form.
//
//  M	        - Input matrix
//  Diagonal	- Output vector of diagonal elements
//  OffDiagonal - Output vector of off diagonal elements
//

template <class T, class M> void TriDiagonal(M& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();
  T       f,g,hh;

  typename M        ::Subscriptor a (A);
  typename Vector<T>::Subscriptor d (   Diagonal);
  typename Vector<T>::Subscriptor od(OffDiagonal);

  for (index_t i=n;i>=2;i--)
  {
    index_t l=i-1;
    T h=0.0, scale=0.0;
    if (l > 1)
    {
      for (index_t k=1;k<=l;k++) scale += fabs(a(i,k));
      if (scale == 0.0)
      {	od(i)=a(i,l);}
      else
      {
	for (index_t k=1;k<=l;k++)
	{
	  a(i,k) /= scale;
	  h += a(i,k)*a(i,k);
	}
	f=a(i,l);
	g = f>0 ? -sqrt(h) : sqrt(h);
	od(i)=scale*g;
	h -= f*g;
	a(i,l)=f-g;
	f=0.0;
	for (index_t j=1;j<=l;j++)
	{
	  a(j,i)=a(i,j)/h;
	  g=0.0;
	  for (index_t k=1  ;k<=j;k++) g += a(j,k)*a(i,k);
	  for (index_t k=j+1;k<=l;k++) g += a(k,j)*a(i,k);
	  od(j)=g/h;
	  f += od(j)*a(i,j);
	}
	hh=f/(h+h);
	for (index_t j=1;j<=l;j++)
	{
	  f=a(i,j);
	  od(j)=g=od(j)-hh*f;
	  for (index_t k=1;k<=j;k++) a(j,k) -= (f*od(k)+g*a(i,k));
	}
      }
    }
    else od(i)=a(i,l);
    d(i)=h;
  }
  d(1)=0.0;
  od(1)=0.0;
  for (index_t i=1;i<=n;i++)
  {
    index_t l=i-1;
    if (d(i))
    {
      for (index_t j=1;j<=l;j++)
      {
	g=0.0;
	for (index_t k=1;k<=l;k++) g += a(i,k)*a(k,j);
	for (index_t k=1;k<=l;k++) a(k,j) -= g*a(k,i);
      }
    }
    d(i)=a(i,i);
    a(i,i)=1.0;
    for (index_t j=1;j<=l;j++) a(j,i)=a(i,j)=0.0;
  }
}

using Type=double;
template void Cholsky<Type>(Matrix<Type>& A);
template void EigenSort(Matrix<Type>&, Vector<Type>&);
template  Matrix<Type> InvertSquare   (const  Matrix<Type>&);
template  Matrix<Type> InvertSymmetric(const  Matrix<Type>&);
template SMatrix<Type> InvertSymmetric(const SMatrix<Type>&);
template void InvertTriangular(Matrix<Type>&);
template void InvertTriangular(SMatrix<Type>&);
template void LUBackSub(const Matrix<Type>&, Vector<Type>&,const std::vector<index_t>&);
template void LUBackSub(const Matrix<Type>&, Matrix<Type>&, const std::vector<index_t>&);
template bool LUDecomp(Matrix<Type>&,std::vector<index_t>&,Type&);
template void QLDecomp(Matrix<Type>&,Vector<Type>&,Vector<Type>&);
template void SVBackSub(const Matrix<Type>&,const Vector<Type>&,const Matrix<Type>&,const Vector<Type>&, Vector<Type>&);
template void SVDecomp<Type,Matrix<Type>>(Matrix<Type>&,Vector<Type>&,Matrix<Type>&);
template void TriDiagonal(Matrix<Type>&,Vector<Type>&,Vector<Type>&);




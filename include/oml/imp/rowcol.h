// File: RowCol.h  Interface base class for RowCol objects.
#ifndef _RowCol_h_
#define _RowCol_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/vecindex.h"

//---------------------------------------------------------------------
//
//  Matrix row proxy class.
//
template <class T, class A, Store M, Data D> class MatrixRow
    : public Indexable<T,MatrixRow<T,A,M,D>,Full,Abstract,VectorShape>
{
public:
    typedef Indexable<T,MatrixRow<T,A,M,D>,Full,Abstract,VectorShape> IndexableT;
    typedef Ref<T,IndexableT,VectorShape> RefT;

    MatrixRow(Indexable<T,A,M,D,MatrixShape>& m,index_t row)
        : itsMatrix(static_cast<A&>(m))
        , itsRow   (row)
    {};
    MatrixRow(const Indexable<T,A,M,D,MatrixShape>& m,index_t row)
        : itsMatrix(static_cast<A&>(const_cast<Indexable<T,A,M,D,MatrixShape>&>(m)))
        , itsRow   (row)
    {};
    MatrixRow(const MatrixRow& r)
        : itsMatrix(r.itsMatrix)
        , itsRow   (r.itsRow)
    {};
    ~MatrixRow() {};
    MatrixRow& operator=(const MatrixRow& r)
    {
        VectorAssign(*this,r);
        return *this;
    }
    MatrixRow& operator=(T t)
    {
        VectorAssign(*this,t);
        return *this;
    }
    template <class A1,Store M1, Data D1> MatrixRow& operator=(const Indexable<T,A1,M1,D1,VectorShape>& v)
    {
        VectorAssign(*this,v);
        return *this;
    }

    index_t   size  () const
    {
        return itsMatrix.GetLimits().Col.size();
    }
    VecLimits GetLimits() const
    {
        return itsMatrix.GetLimits().Col;
    };

    T  operator()(index_t col) const
    {
        return itsMatrix(itsRow,col);
    }
    T& operator()(index_t col)
    {
        return itsMatrix(itsRow,col);
    }

    MatrixRow& operator*=(const T& t)
    {
        for (index_t c:itsMatrix.cols())
            (*this)(c)*=t;
        return *this;
    }

    class Subscriptor
    {
    public:
        Subscriptor(Indexable<T,MatrixRow<T,A,M,D>,Full,Abstract,VectorShape>& a)
            : itsRow(static_cast<MatrixRow<T,A,M,D>&>(a).itsRow)
            , itsS(static_cast<MatrixRow<T,A,M,D>&>(a).itsMatrix)
        {};

        T& operator()(index_t i)
        {
            return itsS(itsRow,i);    //Let matrix do bounds check.
        }

    private:
        const index_t itsRow;
        typename A::Subscriptor itsS;
    };

private:
    friend class Indexable<T,MatrixRow,Full,Abstract,VectorShape>;
    friend class Subscriptor;

    A&    itsMatrix;
    const index_t itsRow;
};



//---------------------------------------------------------------------
//
//  Matrix column proxy class.
//
template <class T, class A, Store M, Data D> class MatrixColumn
    : public Indexable<T,MatrixColumn<T,A,M,D>,Full,Abstract,VectorShape>
{
public:
    typedef Indexable<T,MatrixColumn<T,A,M,D>,Full,Abstract,VectorShape> IndexableT;
    typedef Ref<T,IndexableT,VectorShape> RefT;
    MatrixColumn(Indexable<T,A,M,D,MatrixShape>& m,index_t col)
        : itsMatrix(static_cast<A&>(m))
        , itsColumn(col)
    {};
    MatrixColumn(const Indexable<T,A,M,D,MatrixShape>& m,index_t col)
        : itsMatrix(static_cast<A&>(const_cast<Indexable<T,A,M,D,MatrixShape>&>(m)))
        , itsColumn(col)
    {};
    MatrixColumn(const MatrixColumn& r)
        : itsMatrix(r.itsMatrix)
        , itsColumn(r.itsColumn)
    {};
    ~MatrixColumn() {};
    MatrixColumn& operator=(const MatrixColumn& r)
    {
        VectorAssign(*this,r);
        return *this;
    }
    MatrixColumn& operator=(T t)
    {
        VectorAssign(*this,t);
        return *this;
    }
    template <class A1,Store M1, Data D1> MatrixColumn& operator=(const Indexable<T,A1,M1,D1,VectorShape>& v)
    {
        VectorAssign(*this,v);
        return *this;
    }

    index_t   size  () const
    {
        return itsMatrix.GetLimits().Row.size();
    }
    VecLimits GetLimits() const
    {
        return itsMatrix.GetLimits().Row;
    };

    T  operator()(index_t row) const
    {
        return itsMatrix(row,itsColumn);
    }
    T& operator()(index_t row)
    {
        return itsMatrix(row,itsColumn);
    }

    MatrixColumn& operator*=(const T& t)
    {
        for (index_t r:itsMatrix.rows())
            (*this)(r)*=t;
        return *this;
    }

    class Subscriptor
    {
    public:
        Subscriptor(Indexable<T,MatrixColumn<T,A,M,D>,Full,Abstract,VectorShape>& a)
            : itsColumn(static_cast<MatrixColumn<T,A,M,D>&>(a).itsColumn)
            , itsS(static_cast<MatrixColumn<T,A,M,D>&>(a).itsMatrix)
        {};

        T& operator()(index_t i)
        {
            return itsS(i,itsColumn);    //Let matrix do bounds check.
        }

    private:
        const index_t itsColumn;
        typename A::Subscriptor itsS;
    };

private:
    friend class Indexable<T,MatrixColumn,Full,Abstract,VectorShape>;
    friend class Subscriptor;

    A&    itsMatrix;
    const index_t itsColumn;
};


//---------------------------------------------------------------------
//
//  Matrix diagonal proxy class.
//
template <class T, class A, Store M, Data D> class MatrixDiagonal
    : public Indexable<T,MatrixDiagonal<T,A,M,D>,Full,Abstract,VectorShape>
{
public:
    typedef Indexable<T,MatrixDiagonal<T,A,M,D>,Full,Abstract,VectorShape> IndexableT;
     typedef Ref<T,IndexableT,VectorShape> RefT;

    MatrixDiagonal(Indexable<T,A,M,D,MatrixShape>& m)
        : itsMatrix(static_cast<A&>(m))
    {
//        assert(m.GetLimits().Row==m.GetLimits().Col);
    };
    MatrixDiagonal(const Indexable<T,A,M,D,MatrixShape>& m)
        : itsMatrix(static_cast<A&>(const_cast<Indexable<T,A,M,D,MatrixShape>&>(m)))
    {
//        assert(m.GetLimits().Row==m.GetLimits().Col);
    };
    MatrixDiagonal(const MatrixDiagonal& m)
        : itsMatrix(m.itsMatrix)
    {};
    ~MatrixDiagonal() {};
    MatrixDiagonal& operator=(const MatrixDiagonal& r)
    {
        VectorAssign(*this,r);
        return *this;
    }
    MatrixDiagonal& operator=(T t)
    {
        VectorAssign(*this,t);
        return *this;
    }
    template <class A1,Store M1, Data D1> MatrixDiagonal& operator=(const Indexable<T,A1,M1,D1,VectorShape>& v)
    {
        VectorAssign(*this,v);
        return *this;
    }

    index_t   size  () const
    {
        return Min(itsMatrix.GetLimits().Row.size(),itsMatrix.GetLimits().Col.size());
    }
    VecLimits GetLimits() const
    {
        return itsMatrix.GetLimits().Row.size() >= itsMatrix.GetLimits().Col.size()
            ? itsMatrix.GetLimits().Col
            : itsMatrix.GetLimits().Row;
    };

    T  operator()(index_t i) const
    {
        return itsMatrix(i,i);
    }
    T& operator()(index_t i)
    {
        return itsMatrix(i,i);
    }

    class Subscriptor
    {
    public:
        Subscriptor(Indexable<T,MatrixDiagonal<T,A,M,D>,Full,Abstract,VectorShape>& a)
            : itsS(static_cast<MatrixDiagonal<T,A,M,D>&>(a).itsMatrix)
        {};

        T& operator()(index_t i)
        {
            return itsS(i,i);    //Let matrix do bounds check.
        }

    private:
        typename A::Subscriptor itsS;
    };

private:
    friend class Indexable<T,MatrixDiagonal,Full,Abstract,VectorShape>;
    friend class Subscriptor;

    A&    itsMatrix;
};


#endif //_RowCol_H_

/**
 * _Matrix class
 *      by igor
 *
 * Requires the template parameter to support all of these:
 *  - constructor T( int )
 *  - constructor T( double )
 *  - operator=( T, const T& )
 *  - bool operator<( const T&, const T& )
 *  - T operator-( const T& )
 *  - T operator-=( const T&, const T& )
 *  - T operator*( const T&, const T& )
 *  - T operator/( const T&, const T& )
 *  - operator>>( istream, T )
 *  - operator<<( ostream, const T& )
 *
 * Has
 *  - LU-decomposition: _Matrix::lu()
 *  - Determinant:      _Matrix::det()
 *
 * This file is part of my library of algorithms found here:
 *      http://shygypsy.com/tools/
 * LICENSE:
 *      http://shygypsy.com/tools/LICENSE.html
 * Copyright (c) 2003
 * Contact author:
 *      abednego at gmail.com
 **/
#ifndef _Matrix_CPP
#define _Matrix_CPP

#include <iostream>
#include <string.h>
using namespace std;
template< class T >
class _Matrix
{
    private:
        T *data;
        T *LU;
        int *P, Psgn;

        // What is considered a zero
        static const T EPS = T( 0.0000001 );

    public:
        int m, n;
	_Matrix() { m = n = 0; data = LU = NULL; P = NULL; }
        _Matrix( int m, int n );
        _Matrix( int m, int n, T deflt );
        _Matrix( const _Matrix< T > &mtx );
	_Matrix( const T mtx[], int nRow, int nCol );
	_Matrix( const vector<T> mtx, int nRow, int nCol);
        const _Matrix< T > &operator=( const _Matrix< T > &mtx );
        ~_Matrix();

        T *operator[]( int i );
        inline T &at( T *v, int i, int j ) { return v[i * n + j]; }
        inline T &at( int i, int j ) { return data[i * n + j]; }
        inline T myabs( const T &t ) { return( t < T( 0 ) ? -t : t ); }


	static void add(_Matrix<T>&A, _Matrix<T> &B, _Matrix<T> &result);
	static void subtract(_Matrix<T>&A, _Matrix<T> &B, _Matrix<T> &result);
	static void multiply(_Matrix<T>&A, _Matrix<T> &B, _Matrix<T> &result);
	static void transpose(_Matrix<T>&A, _Matrix<T> &result);
	static void invUpper(_Matrix<T> &upMat, _Matrix<T> &result);
	static void invLower(_Matrix<T> &lowMat, _Matrix<T> &result);
	static void inverse(_Matrix<T> &A, _Matrix<T> & result);

        /**
         * The PLU-decomposition is stored in two companion arrays
         * P and LU. They are allocated by lu() and become obsolete
         * whenever the original matrix is modified. It is the user's
         * responsibility to call lu() again when that happens. The
         * () operator is used to access the LU matrix, just like []
         * is used for the matrix data. p() accesses the permutation
         * vector P.
         **/
        void lu();
        T &operator()( int i, int j );
        int &p( int i );

    // I/O Friends
    //friend istream &operator>><>( istream &in, _Matrix< T > &mtx );
    //friend ostream &operator<<<>( ostream &out, const _Matrix< T > &mtx );

    // Other friends
    //friend T det<>( _Matrix< T > &mtx );
    T det();

    // DEBUG
    void printLU()
    {
        for( int i = 0; i < m; i++ )
        {
            for( int j = 0; j < n; j++ )
            {
                if( j ) cout << " ";
                cout << at( LU, i, j );
            }
            cout << endl;
        }
    }

    void printP()
    {
	if (P==NULL) return;
	for( int i = 0; i < m; i++)
	   cout << P[i] << " " ;
	cout << endl;
    }
};

template< class T > istream &operator>>( istream &in, _Matrix< T > &mtx );
template< class T > ostream &operator<<( ostream &out, const _Matrix< T > &mtx );

//------------------ Implementation ---------------------//
template< class T >
_Matrix< T >::_Matrix( int m, int n )
{
    this->m = m;
    this->n = n;
    data = new T[m * n];
    LU = NULL;
    P = NULL;
}

template< class T >
_Matrix< T >::_Matrix( int m, int n, T deflt )
{
    this->m = m;
    this->n = n;
    data = new T[m * n];
    LU = NULL;
    P = NULL;

    for( int i = 0; i < m * n; i++ )
        data[i] = deflt;
}

template< class T > _Matrix< T >::_Matrix( const _Matrix< T > &mtx )
{
    this->m = mtx.m;
    this->n = mtx.n;
    data = new T[m * n];
    LU = NULL;
    P = NULL;


    for( int i = 0; i < m * n; i++ )
        data[i] = mtx.data[i];
}

template< class T > _Matrix< T >::_Matrix( const T mtx[], int nRow, int nCol )
{
    this->m = nRow;
    this->n = nCol;
    data = new T[m * n];
    LU = NULL;
    P = NULL;


    for( int i = 0; i < m * n; i++ )
        data[i] = mtx[i];
}

template< class T > _Matrix< T >::_Matrix( const vector<T> mtx, int nRow, int nCol)
{
    this->m = nRow;
    this->n = nCol;
    data = new T[m * n];
    LU = NULL;
    P = NULL;


    for( int i = 0; i < m * n; i++ )
        data[i] = mtx[i];
}

template< class T > void _Matrix< T >::add(_Matrix<T> &A, _Matrix<T> &B, _Matrix<T> &result)
{
	for(int i=0; i<A.m; i++)
	{
		for(int j=0; j<A.n; j++) result.at(i, j) = A.at(i, j) + B.at(i,  j);
	}
}

template< class T > void _Matrix< T >::subtract(_Matrix<T> &A, _Matrix<T> &B, _Matrix<T> &result)
{
	for(int i=0; i<A.m; i++)
	{
		for(int j=0; j<A.n; j++) result.at(i, j) = A.at(i, j) - B.at(i,  j);
	}
}
	
template< class T> void _Matrix<T>::multiply(_Matrix<T>&A, _Matrix<T> &B, _Matrix<T> &result)
{
	// m = row, n = col
	if (A.n != B.m) 
	{
		cout << "Dimension mismatched" << endl;
		return;
	}
	int C = A.n;
	for (int i = 0; i < A.m; i++)
		for (int j = 0; j < B.n; j++) 
		{
			T res = 0;
			for (int c = 0; c < C; c++)
				res += A.at(i,c) * B.at(c,j);
			result.at(i,j) = res;
		}
}
template <class T> void _Matrix<T>::transpose(_Matrix<T>&A, _Matrix<T> &result)
{
	if (A.n != result.m || A.m != result.n)
	{
		cout << "Dimension mismatched" << endl;
		return;
	}
	for (int i =0; i < A.m; i++)
		for (int j =0; j < A.n; j++)
			result.at(j,i) = A.at(i,j);
}
template< class T > void  _Matrix< T >::invUpper(_Matrix<T> &upMat, _Matrix<T> &result)
{
	// assuming square matrix
	int n = upMat.m;

	// we solve column by column, from right to elft
	for(int i=n-1; i>=0; i--)
	{		
		result.at(i, i) = 1 / upMat.at(i, i);
		// from bottom to top
		for(int j=i-1; j>=0; j--)
		{
			T sum = 0;
			for(int k=n-1; k>j; k--) sum-=upMat.at(j, k)*result.at(k, i);
			result.at(j, i) = sum / upMat.at(j, j);
		}		
	}
}

template< class T > void  _Matrix< T >::invLower(_Matrix<T> &lowMat, _Matrix<T> &result)
{
        // assuming square matrix
        int n = lowMat.m;

        // we solve column by column, from right to elft
        for(int i=0; i<n; i++)
        {
                result.at(i, i) = 1 / lowMat.at(i, i);
                // from bottom to top
                for(int j=i+1; j<n; j++)
                {
                        T sum = 0;
                        for(int k=0; k<j; k++) sum-=lowMat.at(j, k)*result.at(k, i);
                        result.at(j, i) = sum / lowMat.at(j, j);
                }
        }
}


template< class T > void _Matrix<T>::inverse(_Matrix<T> &mtx, _Matrix<T> & result)
{
    	if( mtx.m != mtx.n ) 
	{	
		throw "Not a square matrix";
		return;
	}

	_Matrix<T> pMat(mtx.m, mtx.n, 0);
    	_Matrix<T> lowMat(mtx.m, mtx.n, 0);
	_Matrix<T> upMat(mtx.m, mtx.n, 0);
	mtx.lu();
	// Copy from LU to lowMat and upMat
 	for (int i = (mtx.m-1); i >=0; i--)
	{
		for (int j = 0; j <=i; j++)
			if (i == j ) 
				lowMat.at(i, j) = 1;
			else
				lowMat.at(i, j) = mtx(i, j);
	}

	for (int i = 0; i < mtx.m; i++)
	{
		for (int j = i; j < mtx.n; j++)
			upMat.at(i, j) = mtx(i,j);
	}
	// Copy from P to pMat
    	for (int i = 0; i < mtx.m; i++)
		pMat.at(i,mtx.p(i)) = 1;
	
	// Inverse _Matrix
	_Matrix<T> tempMat1(mtx.m, mtx.n);
	_Matrix<T> tempMat2(mtx.m, mtx.n);
	
	_Matrix<T>::invUpper(upMat, tempMat1);
	_Matrix<T>::invLower(lowMat, tempMat2);
	_Matrix<T>::multiply(tempMat1, tempMat2, result);

	tempMat1 = result;
//	_Matrix<T>::transpose(pMat, tempMat2);
	_Matrix<T>::multiply(tempMat1, pMat, result);
}

template< class T >
const _Matrix< T > &_Matrix< T >::operator=( const _Matrix< T > &mtx )
{
    if( &mtx != this )
    {
        delete [] data;
        m = mtx.m;
        n = mtx.n;
        data = new T[m * n];
        for( int i = 0; i < m * n; i++ )
            data[i] = mtx.data[i];
    }
    return *this;
}

template< class T >
_Matrix< T >::~_Matrix()
{
    delete [] data;
}

template< class T >
T *_Matrix< T >::operator[]( int i )
{
    return data + i * n;
}

template< class T >
void _Matrix< T >::lu()
{
    if( LU ) delete [] LU;
    if( P ) delete [] P;
    LU = new T[m * n];
    P = new int[m];
    memcpy( LU, data, m * n * sizeof( T ) );
    for( int i = 0; i < m; i++ ) P[i] = i;
    Psgn = 1;

    for( int r = 0, c = 0; r < m && c < n; r++, c++ )        // For each row
    {
        // Find largest pivot in this column
        int pr = r;
        for( int i = r + 1; i < m; i++ )
            if( myabs( at( LU, i, c ) ) > myabs( at( LU, pr, c ) ) )
                pr = i;
        if( myabs( at( LU, pr, c ) ) <= EPS )
        {
            // Singular matrix; skip column
            r--;
            continue;
        }

        if( pr != r )
        {
            // Swap rows r and pr
            P[r] ^= P[pr] ^= P[r] ^= P[pr];
            Psgn = -Psgn;
            for( int i = 0; i < n; i++ )
            {
                T tmp = at( LU, r, i );
                at( LU, r, i ) = at( LU, pr, i );
                at( LU, pr, i ) = tmp;
            }
        }

        // Subtract row r from rows below it
        for( int s = r + 1; s < m; s++ )
        {
            at( LU, s, c ) = at( LU, s, c ) / at( LU, r, c );
            for( int d = c + 1; d < n; d++ )
                at( LU, s, d ) -= at( LU, s, c ) * at( LU, r, d );
        }
    }
}

template< class T >
T &_Matrix< T >::operator()( int i, int j )
{
    return LU[i * n + j];
}

template< class T >
int &_Matrix< T >::p( int i )
{
    return P[i];
}

template< class T >
istream &operator>>( istream &in, _Matrix< T > &mtx )
{
    for( int i = 0; i < mtx.m * mtx.n; i++ )
        in >> mtx.data[i];
    return in;
}

template< class T >
ostream &operator<<( ostream &out, _Matrix< T > &mtx )
{
    for( int i = 0; i < mtx.m; i++ )
    {
        for( int j = 0; j < mtx.n; j++ )
        {
            if( j ) out << " ";
            out << mtx.at( i, j );
        }
        out << endl;
    }
    return out;
}

template< class T >
T _Matrix<T>::det()
{
    if( m != n ) throw "Not a square matrix";
    lu();
    T ans = 1;
    for( int i = 0; i < n; i++ )
        ans *= at( LU, i, i );
    return( Psgn > 0 ? ans : -ans );
}


//------------------- DEBUG and Testing ------------------//
/* 
int main()
{
    int x, y;
    cin >> x >> y;
    _Matrix< double > m( x, y );
    cin >> m;
    cout << "det(M) = " << det( m ) << endl;
    return 0;
}
*/
#endif

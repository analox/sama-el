/*! RBFN-PTO.cpp
* ======================================================================
*
*  File    :  RBFN-PTO.cpp
*  Created :  01.01.1995
*
*  Copyright (c) 1995,1998 Martin Kreutz
*
*
*   \par Copyright (c) 1999-2003:
*        Institut f&uuml;r Neuroinformatik<BR>
*        Ruhr-Universit&auml;t Bochum<BR>
*        D-44780 Bochum, Germany<BR>
*        Phone: +49-234-32-25558<BR>
*        Fax:   +49-234-32-14209<BR>
*        eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
*        www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
*        <BR>
* 
*    \par Project:
*        Mixture
*  
*    \par Language and Compiler:
*        C++, egcs (Linux)
*  
*    \par File and Revision:
*        $RCSfile: RBFN-PTO.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: RBFN-PTO.cpp,v $
*        Revision 2.1  2004/06/18 00:50:11  saviapbe
*        The standard file header was adapted for doxygen.
*
*        Revision 2.0  2003/11/28 16:23:12  shark-admin
*        Revision tag reset to revision tag 2.x
*  
*        Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
*        INI Administration
*  <BR>
* 
*  Last update: 19.10.1998
*
* ----------------------------------------------------------------------
*
*  This file is part of the EALib. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  This library is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this library; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
* ======================================================================
*/

#include "EALib/sqr.h"
#include "Array/ArrayOp.h"
#include "Mixture/RBFN-PTO.h"

const double MIN_VAL = 1e-100;

//===========================================================================

RBFN_PTO::RBFN_PTO( unsigned nInputs,
		    unsigned nOutputs,
		    unsigned nCenters,
		    bool     trainkap )
  : RBFN( nInputs, nOutputs, nCenters ),
    kappa     ( nCenters ),
    dkappa    ( nCenters ),
    trainkappa( trainkap )
{
    kappa = 1;
}

//===========================================================================

void RBFN_PTO::copyFrom( const array< double >& w )
{
    SIZE_CHECK( w.nelem( ) == b.nelem( ) +
		              RBFN::A.nelem( ) +
		              m.nelem( ) +
		              v.nelem( ) +
		              kappa.nelem( ) )

    unsigned i, k;

    for( k = i = 0; i < b.nelem( ); ++i, ++k )
        b.elem( i ) = w( k );
    for( i = 0; i < RBFN::A.nelem( ); ++i, ++k )
        RBFN::A.elem( i ) = w( k );
    for( i = 0; i < m.nelem( ); ++i, ++k )
        m.elem( i ) = w( k );
    for( i = 0; i < v.nelem( ); ++i, ++k )
        v.elem( i ) = sqr( w( k ) );
    for( i = 0; i < kappa.nelem( ); ++i, ++k )
        kappa.elem( i ) = w( k );
}

void RBFN_PTO::copyTo( array< double >& w ) const
{
    unsigned i, k;

    w.resize( b.nelem( ) +
	      RBFN::A.nelem( ) +
	      m.nelem( ) +
	      v.nelem( ) +
	      kappa.nelem( ) );

    for( k = i = 0; i < b.nelem( ); ++i, ++k )
        w( k ) = b.elem( i );
    for( i = 0; i < RBFN::A.nelem( ); ++i, ++k )
        w( k ) = RBFN::A.elem( i );
    for( i = 0; i < m.nelem( ); ++i, ++k )
        w( k ) = m.elem( i );
    for( i = 0; i < v.nelem( ); ++i, ++k )
        w( k ) = sqrt( v.elem( i ) );
    for( i = 0; i < kappa.nelem( ); ++i, ++k )
        w( k ) = kappa.elem( i );
}

void RBFN_PTO::copyGrad( array< double >& dw ) const
{
    unsigned i, k;

    dw.resize( dB.nelem( ) +
	       dA.nelem( ) +
	       dm.nelem( ) +
	       ds.nelem( ) +
	       dkappa.nelem( ) );

    for( k = i = 0; i < dB.nelem( ); ++i, ++k )
        dw( k ) = dB.elem( i );
    for( i = 0; i < dA.nelem( ); ++i, ++k )
        dw( k ) = dA.elem( i );
    for( i = 0; i < dm.nelem( ); ++i, ++k )
        dw( k ) = dm.elem( i );
    for( i = 0; i < ds.nelem( ); ++i, ++k )
        dw( k ) = ds.elem( i );
    for( i = 0; i < dkappa.nelem( ); ++i, ++k )
        dw( k ) = dkappa.elem( i );
}

//===========================================================================

void RBFN_PTO::copyFrom( const MixtureOfGaussians& gm )
{
    unsigned i, j, k;
    double   varg;

    for( i = 0; i < size( ); ++i ) {
	for( varg = 1, j = 0; j < dim( ); ++j ) {
	    m( i, j ) = gm.m( i, j );
	    varg *= ( v( i, j ) = std::max( gm.v( i, j ), MIN_VAL ) );
	}

	for( k = 0; j < gm.dim( ); ++j, ++k )
	    RBFN::A( k, i ) = gm.m( i, j );

	kappa( i ) = gm.a( i ) / sqrt( varg );
    }

    b = 0;

    clear( );
}

//===========================================================================

void RBFN_PTO::firstLayer( const array< double >& x,
			   array< double >& ex ) const
{
    SIZE_CHECK( x.ndim( ) == 1 && x.dim( 0 ) == dim( ) )

    unsigned i, j;
    double   marg, norm;

    ex.resize( size( ) );

    for( norm = 0, i = 0; i < size( ); i++ ) {
        for( marg = 0, j = dim( ); j--; )
	    marg -= sqr( x( j ) - m( i, j ) ) / std::max( v( i, j ), MIN_VAL );
	norm += ( ex( i ) = kappa( i ) * exp( marg / 2 ) );
    }

    ex /= norm;
}

//===========================================================================

void RBFN_PTO::gradientClear( )
{
    dkappa.resize( kappa );
    dkappa = 0.;
    RBFN::gradientClear( );
}

//===========================================================================

//
// nur fuer 1-dim. Ausgang !!!
//
void RBFN_PTO::gradientMSE( )
{
    unsigned i, k, l;
    double s, t;
    array< double > e, Ae;
    array< double > tdB, tdA, tdk;
    array_reference< double > x, y;

    for( l = 0; ; ++l ) {
        if( _x->ndim( ) == 2 ) {
	    x.copyReference( ( *_x )[ l ] );
	    y.copyReference( ( *_y )[ l ] );
	} else {
	    x.copyReference( *_x );
	    y.copyReference( *_y );
	}

	firstLayer( x, e );

	Ae  = inner_product( RBFN::A, e );
	dB += ( tdB = ( Ae + b - y ) * 2. );
	dA += ( tdA = outer_product( tdB, e ) );
	tdk.resize( size( ) );

	for( k = 0; k < size( ); ++k ) {
	    dkappa( k ) += ( tdk( k ) = dA( 0, k ) * ( RBFN::A( 0, k ) - Ae( 0 ) ) ) / kappa( k );
	    for( i = 0; i < dim( ); ++i ) {
	        s  = sqrt( v( k, i ) );
		t  = x( i ) - m( k, i );

		dm( k, i ) += tdk( k ) * t / ( s*s );
		ds( k, i ) += tdk( k ) * ( t*t ) / ( s*s*s );
	    }
	}

	err += sumOfSqr( tdB / 2. );

	if( _x->ndim( ) == 1 || l+1 >= _x->dim( 0 ) )
	    break;
    }
}

//===========================================================================

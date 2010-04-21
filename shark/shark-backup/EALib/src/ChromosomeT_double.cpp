/*! ChromosomeT_double.cpp
* ======================================================================
*
*  File    :  ChromosomeT_double.cpp
*  Created :  1995-01-01
*
*  Copyright (c) 1995-2000 Martin Kreutz
*  Copyright (c) 1999-2000 Bernhard Sendhoff

*  Institut fuer Neuroinformatik
*  Ruhr-Universitaet Bochum
*  44780 Bochum, Germany<BR>
*  Phone: +49-234-32-25558<BR>
*  Fax:   +49-234-32-14209<BR>
*  eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
*  www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
*  
*    \par Project:
*        EALib
*  
*    \par Language and Compiler:
*        C++, egcs (Linux)
*  
*    \par File and Revision:
*        $RCSfile: ChromosomeT_double.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: ChromosomeT_double.cpp,v $
*        Revision 2.7  2005/02/03 18:08:22  christian_igel
*        *** empty log message ***
*
*        Revision 2.6  2005/02/03 18:03:03  christian_igel
*        real-coded GA operators updated
*
*        Revision 2.5  2004/12/15 15:31:00  shark-admin
*        #ifdef for win32 C1001 error deactivated
*
*        Revision 2.4  2004/12/10 18:33:11  christian_igel
*        Real-valued GA operators.
*
*        Revision 2.3  2004/06/17 23:08:59  saviapbe
*        An error in the INT's URL was corrected.
*
*        Revision 2.2  2004/06/17 22:53:19  saviapbe
*        The standard file header was for doxygen adapted.
*
*        Revision 2.1  2004/05/26 19:50:15  saviapbe
*  
*  
*        A bug with the multiplication of rotations matrices in the method "ChromosomeT<double>::mutateRotate" was eliminated.
*  
*        Revision 2.0  2003/11/28 16:23:09  shark-admin
*        Revision tag reset to revision tag 2.x
*  
*        Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
*        INI Administration
*  <BR>
*  
*  
*  Last update: 18th, Feb, 2001 by Tatsuya Okabe (HONDA R&D EUROPE)
*       tatsuya.okabe@hre-ftr.f.rd.honda.co.jp
*       Function mutateRotate( ChromosomeT<double>, double, double,
*                double, int, double ) sign_of_sigma = 0 -> -1
*
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
#include "LinAlg/linalg.h"
//#include "LinAlg/arraylinalg.h"
#include "EALib/ChromosomeT.h"

//===========================================================================

ChromosomeT< double >::~ChromosomeT( )
{
}

//===========================================================================

void ChromosomeT< double >::initializeDerandom( double minVal, double maxVal )
{
    SIZE_CHECK( size( ) >= 5 )

    unsigned n = ( size( ) - 2 ) / 3;
    std::vector< double >::iterator s  = begin( );
    std::vector< double >::iterator d  = s + n;
    std::vector< double >::iterator r  = d + n;
    std::vector< double >::iterator dr = r + n;
    std::vector< double >::iterator sr = dr+ 1;

    for( unsigned i = 0; i < n; i++ ) {
        //
        // kumulierte Zufallsschritte
        //
        *( s + i ) = 0;

	//
	// Schrittweite fuer Mutation der Objektvariablen
	//
	*( d + i ) = Rng::uni( minVal, maxVal );;

	//
	// zu lernende Vorzugsrichtung
	//
	*( r + i ) = 0;
    }

    //
    // Schrittweite fuer Mutation der Vorzugsrichtung
    //
    *dr = Rng::uni( minVal, maxVal );

    //
    // kumulierte Schritte (Richtung)
    //
    *sr = 0;
}

//===========================================================================

void ChromosomeT< double >::dumpDerandom( std::ostream& os ) const
{
    SIZE_CHECK( size( ) >= 5 )

    unsigned i, n = ( size( ) - 2 ) / 3;
    std::vector< double >::const_iterator s  = begin( );
    std::vector< double >::const_iterator d  = s + n;
    std::vector< double >::const_iterator r  = d + n;
    std::vector< double >::const_iterator dr = r + n;
    std::vector< double >::const_iterator sr = dr+ 1;

    os << "s =";
    for( i = 0; i < n; i++ )
        os << '\t' << *( s + i );
    os << "\nd =";
    for( i = 0; i < n; i++ )
        os << '\t' << *( d + i );
    os << "\nr =";
    for( i = 0; i < n; i++ )
        os << '\t' << *( r + i );
    os << "\ndr =\t" << *dr;
    os << "\nsr =\t" << *sr << std::endl;
}

//===========================================================================

inline unsigned long pow2( unsigned n )
{
    return 1UL << n;
}

void ChromosomeT< double >::decodeBinary( const Chromosome& chrom,
				          const Interval&   range,
				          unsigned          nbits,
				          bool              useGray )
{
    SIZE_CHECK( chrom.size( ) % nbits == 0 )

    const std::vector< bool >& src = dynamic_cast< const std::vector< bool >& >( chrom );

    double stepSize = ( pow2( nbits ) - 1 ) / range.width();
    unsigned i, j, k;
    unsigned long l, m;

    resize( src.size( ) / nbits );


    for( j = size( ), k = src.size( ); j--; ) {
	if( useGray )
	    for( l = m = 0, i = nbits; i--; )
                l = ( l << 1 ) | ( m ^= ( src[ --k ] ? 1 : 0 ) );
	else
	    for( l = 0, i = nbits; i--; )
                l = ( l << 1 ) | ( src[ --k ] ? 1 : 0 );

	( *this )[ j ] = l / stepSize + range.lowerBound();
    }
}

//===========================================================================

void ChromosomeT< double >::accumulate( const std::vector< double >& acc, double c )
{
    SIZE_CHECK( size( ) == acc.size( ) )

    for( unsigned i = size( ); i--; )
        ( *this )[ i ] = ( 1 - c ) * ( *this )[ i ] + c * acc[ i ];
}

//===========================================================================

void ChromosomeT< double >::accumulate( const Chromosome& acc, double c )
{
    accumulate( dynamic_cast< const std::vector< double >& >( acc ), c );
}

//===========================================================================

void ChromosomeT< double >::mutateNormal( double variance )
{
    for( unsigned i = size( ); i--; )
        ( *this )[ i ] += Rng::gauss( 0, variance );
}

//===========================================================================

void ChromosomeT< double >::mutateNormal( const std::vector< double >& variances, bool cycle )
{
    RANGE_CHECK( variances.size( ) <= size( ) )

    for( unsigned i = cycle ? size( ) : variances.size( ); i--; )
        ( *this )[ i ] += Rng::gauss( 0, variances[ i % variances.size( ) ] );
}

//===========================================================================

void ChromosomeT< double >::mutateNormal( const ChromosomeT< double >& variances,
				          bool cycle )
{
    mutateNormal( static_cast< const std::vector< double >& >( variances ), cycle );
}

//===========================================================================

void ChromosomeT< double >::mutateNormal( const Chromosome& variances,
				          bool cycle )
{
    mutateNormal( dynamic_cast< const std::vector< double >& >( variances ), cycle );
}

//===========================================================================

void ChromosomeT< double >::mutateCauchy( double scale )
{
    for( unsigned i = size( ); i--; )
        ( *this )[ i ] += Rng::cauchy( ) * scale;
}

//===========================================================================

void ChromosomeT< double >::mutateCauchy( const std::vector< double >& scale,
				          bool cycle )
{
    RANGE_CHECK( scale.size( ) <= size( ) )

    for( unsigned i = cycle ? size( ) : scale.size( ); i--; )
        ( *this )[ i ] += Rng::cauchy( ) * scale[ i % scale.size( ) ];
}

//===========================================================================

void ChromosomeT< double >::mutateCauchy( const ChromosomeT< double >& scale,
				          bool cycle )
{
    mutateNormal( static_cast< const std::vector< double >& >( scale ), cycle );
}

//===========================================================================

void ChromosomeT< double >::mutateCauchy( const Chromosome& scale,
				          bool cycle )
{
    mutateNormal( dynamic_cast< const std::vector< double >& >( scale ), cycle );
}

//===========================================================================

void ChromosomeT< double >::mutateNormalRotAngles( const Chromosome& sigma,
					           const Chromosome& alpha )
{
    mutateNormalRotAngles( dynamic_cast< const std::vector< double >& >( sigma ),
                           dynamic_cast< const std::vector< double >& >( alpha ) );
}

//===========================================================================

void
ChromosomeT< double >::mutateNormalRotAngles( const std::vector< double >& sigma_sqr,
					      const std::vector< double >& alpha )
{
    RANGE_CHECK( sigma_sqr.size( ) > 0 || size( ) == 0 )

    register unsigned m, i, j, k;

    std::vector< double >& v = *this;
    std::vector< double > z( v.size( ) );

    //
    // create random vector z
    //
    for( m = 0; m < z.size( ) && m < sigma_sqr.size( ); m++ )
        z[ m ] = Rng::gauss( 0, sigma_sqr[ m ] );

    for( ; m < z.size( ); m++ )
        z[ m ] = Rng::gauss( 0, sigma_sqr[ sigma_sqr.size( ) - 1 ] );

    //
    // rotate vector z
    //
    for( k = 0, j = 1; j < m && k < alpha.size( ); j++ )
        for( i = 0; i < m-j && k < alpha.size( ); i++, k++ ) {
	    double sina = sin( alpha[ k ] );
	    double cosa = cos( alpha[ k ] );
	    double t    = cosa * z[ i ] - sina * z[ i + j ];
	    z[ i + j ]  = sina * z[ i ] + cosa * z[ i + j ];
	    z[ i     ]  = t;
	}

    //
    // mutate *this with rotated random vector z
    //
    for( m = v.size( ); m--; )
        v[ m ] += z[ m ];
}

//===========================================================================

void ChromosomeT< double >::mutateLogNormal( double overallVariance,
					     double indivVariance )
{
    double overall = Rng::gauss( 0, overallVariance );

    for( unsigned i = size( ); i--; )
        ( *this )[ i ] *= exp( overall + Rng::gauss( 0, indivVariance ) );
}

//===========================================================================

void ChromosomeT< double >::mutateDerandom( std::vector< double >& v,
				            const DerandomConst& K )
{
    SIZE_CHECK( v.size( ) == 3 * size( ) + 2 )

    unsigned n = size( );
    std::vector< double >::iterator s  = v.begin( );
    std::vector< double >::iterator d  = s + n;
    std::vector< double >::iterator r  = d + n;
    std::vector< double >::iterator dr = r + n;
    std::vector< double >::iterator sr = dr+ 1;

    unsigned i;
    double   z, zr, globalFac, sum;
    std::vector< double > y( *this );

    zr = Rng::gauss( 0, 1 );

    for( sum = 0, i = 0; i < size( ); i++ ) {

	//
	// create random step
	//
	z = Rng::gauss( 0, 1 );

	//
	// mutate object variables
	//
	( *this )[ i ] += *( d + i ) * z + *dr * zr * *( r + i );

	//
	// accumulate step sizes
	//
	*( s + i ) = ( 1-K.C ) * *( s + i ) + K.C * ( K.Cu*z );
	sum += *( s + i ) * *( s + i );
    }

    //
    // step size factor
    //
    globalFac = exp( K.Beta * ( sqrt( sum ) - K.ChiN ) );

    //
    // mutate step sizes
    //
    for( i = 0; i < size( ); i++ )
        *( d + i ) *= globalFac * exp( K.BetaI * ( fabs(*(s+i)) - K.Chi1 ) );

    //
    // accumulate direction
    //
    *sr = std::max( ( 1-K.C ) * *sr + K.C * ( K.Cu*zr ), 0. );

    //
    // mutate step width for direction
    //
    for( sum = 0, i = 0; i < size( ); i++ ) {
	y[ i ] = ( 1-K.Crr ) * *dr * *(r+i) + K.Crr * ( (*this)[i]-y[i] );
	sum += y[ i ] * y[ i ];
    }

    //
    // normalize direction vector
    //
    if( ( sum = sqrt( sum ) ) > 0 )
        for( i = 0; i < size( ); i++ )
	    *( r + i ) = y[ i ] / sum;

    for( sum = 0, i = 0; i < size( ); i++ )
	sum += *( d + i ) * *( d + i );

    sum = sqrt( sum ) / 3;

    *dr *= exp( K.BetaR * ( fabs( *sr ) - K.Chi1 ) );
    if( *dr < sum ) *dr = sum;
}

//===========================================================================

void ChromosomeT< double >::mutateDerandom( Chromosome& chrom,
				            const DerandomConst& K )
{
    mutateDerandom( dynamic_cast< std::vector< double >& >( chrom ), K );
}

//===========================================================================

void
ChromosomeT< double >::recombineIntermediate( const Chromosome& dadChrom,
					      const Chromosome& momChrom )
{
    SIZE_CHECK( dadChrom.size( ) == momChrom.size( ) )

    const std::vector< double >& dad = dynamic_cast< const std::vector< double >& >( dadChrom );
    const std::vector< double >& mom = dynamic_cast< const std::vector< double >& >( momChrom );

    resize( dad.size( ) );

    for( unsigned i = size( ); i--; )
        ( *this )[ i ] = ( dad[ i ] + mom[ i ] ) / 2.;
}

//===========================================================================
//
// generalized intermediate recombination
//
void
ChromosomeT< double >::recombineGenIntermediate( const Chromosome& dadChrom,
					         const Chromosome& momChrom )
{
    SIZE_CHECK( dadChrom.size( ) == momChrom.size( ) )

    const std::vector< double >& dad = dynamic_cast< const std::vector< double >& >( dadChrom );
    const std::vector< double >& mom = dynamic_cast< const std::vector< double >& >( momChrom );

    resize( dad.size( ) );

    for( unsigned i = size( ); i--; )
        ( *this )[ i ] = Rng::uni( dad[ i ], mom[ i ] );
}

//===========================================================================
//
// generalized intermediate recombination
//
void
ChromosomeT< double >::recombineGeomIntermediate( const Chromosome& dadChrom,
					          const Chromosome& momChrom )
{
    SIZE_CHECK( dadChrom.size( ) == momChrom.size( ) )

    const std::vector< double >& dad = dynamic_cast< const std::vector< double >& >( dadChrom );
    const std::vector< double >& mom = dynamic_cast< const std::vector< double >& >( momChrom );

    resize( dad.size( ) );

    for( unsigned i = size( ); i--; )
        ( *this )[ i ] = sqrt( dad[ i ] * mom[ i ] );
}

//===========================================================================

void ChromosomeT< double >::recombineIntermediate( Chromosome& mate )
{
    SIZE_CHECK( size( ) == mate.size( ) )

    std::vector< double >& v = *this;
    std::vector< double >& w = dynamic_cast< std::vector< double >& >( mate );

    for( unsigned i = v.size( ); i--; )
        v[ i ] = w[ i ] = ( v[ i ] + w[ i ] ) / 2.;
}

//===========================================================================
//
// generalized intermediate recombination
//
void ChromosomeT< double >::recombineGenIntermediate( Chromosome& mate )
{
    SIZE_CHECK( size( ) == mate.size( ) )

    std::vector< double >& v = *this;
    std::vector< double >& w = dynamic_cast< std::vector< double >& >( mate );

    for( unsigned i = v.size( ); i--; ) {
	double a = v[ i ];
	double b = w[ i ];
	v[ i ] = Rng::uni( a, b );
	w[ i ] = Rng::uni( a, b );
    }
}

//===========================================================================
//
// generalized intermediate recombination
//
void ChromosomeT< double >::recombineGeomIntermediate( Chromosome& mate )
{
    SIZE_CHECK( size( ) == mate.size( ) )

    std::vector< double >& v = *this;
    std::vector< double >& w = dynamic_cast< std::vector< double >& >( mate );

    for( unsigned i = v.size( ); i--; )
        v[ i ] = w[ i ] = sqrt( v[ i ] * w[ i ] );
}

//===========================================================================

#ifdef OLD_STUFF

//
// cf. B. Sendhoff & M. Kreutz: Analysis of possible genome-dependence of
//     mutation rates in Genetic Algorithms.
//     Evolutionary Computing - selected papers from the 1996 AISB workshop,
//     vol. 1134 of Lecture Notes in Computer Science, T. Fogarty, ed.,
//     pp. 257-269, Springer Verlag, Berlin, September 1996.
//

static double a( long i, double s )
{
    const unsigned Nk = 6;
    static double k[ 2 ][ Nk ] = {
        {  0.5437, -0.3672, -0.0634,  0.6627, -0.5034,  0.0021 },
	{  0.0043,  0.1677, -0.2843,  0.0073,  1.1447, -0.8107 }
    };

    unsigned j;
    double l, logs, sum = 0.0;

    logs = log( 1 + 1 / s ) / Ln2;

    for( j = 0, l = 1.; j < Nk; j++, l *= logs )
	sum += k[ i ][ j ] * l;

    return sum;
}

static inline double q1( double i, double a1, double a2 )
{
    register double pow2mi1 = pow( 2., -i - 1 );

    return a1 * ( 1 - pow2mi1 )
         - a2 * ( i*i - 2 * i + 3 - 6 * pow2mi1 );
}

static double q2( double x, double y )
{
    return 0.5 * x * exp( -y )
         +   2 * x * exp( - 4 * y )
         +   8 * x * exp( -16 * y );
}

static inline double Pd( double i, double a1, double a2 )
{
    double p = ( a1 - a2 * i * i );
    return p > 0 ? p : 0;
}

static inline double Pc( double i, double a1, double a2 )
{
    double sqrt_a1a2 = sqrt( fabs( a1 / a2 ) );

    if( i > sqrt_a1a2 )
	return q1( sqrt_a1a2, a1, a2 ) * exp( -Ln2 * ( i - sqrt_a1a2 ) );
    else
        return q1( i, a1, a2 );
}

static inline double P( double i, double a1, double a2, double step )
{
    if( i == 0 )
        return a1; // Pd( 0, a1, a2 );
    else {
	double pd = Pd( i,      a1, a2 );
	double pc = Pc( i-step, a1, a2 );
	return pd * ( 1 - pc ) + pc * ( 1- pd );
    }
}

void ChromosomeT< double >::stdDevToMutationRates( double stdDev, double rate )
{
    const unsigned NBits = 10;
    register unsigned  i;
    double eta, a1, a2, c1, c2, p1, p2, qq2, id, pow2i, step;

    eta = 0.4577 - 0.5052 * sqr( 1 - stdDev );
    a1  = a( 0, stdDev );
    a2  = a( 1, stdDev );
    c1  = 2 / ( sqrt( 2 * Pi ) * stdDev );
    c2  = sqr( ( 1 + eta ) / ( Sqrt2 * stdDev ) );

    if( stdDev < 1 ) qq2 = q2( c1, c2 );

    std::vector< double >& p = *this;
    step = p.size( ) > NBits
        ? double( NBits ) / p.size( ) : 1.;

    for( i = 0; i < p.size( ); i++ ) {
	id = i * step;
	if( stdDev < 1 ) {
	    if( i == 0 )
	        p[ i ] = c1 * exp( -c2 );
	    else {
		pow2i  = pow( 2., id );
		p1     = c1 * pow2i * exp( -c2 * pow2i * pow2i );
		p2     = qq2 / pow2i;
		p[ i ] = p1 * ( 1 - p2 ) + p2 * ( 1 - p1 );
	    }
	} else {
	    if( i == 0 )
	        p[ i ] = a1;
	    else
	        p[ i ] = P( id, a1, a2, step );
	}
    }

    //
    // scale p[ i ]
    //
    double sum = 0.;
    for( i = 0; i < p.size( ); i++ )
	sum += p[ i ];
    rate /= sum;
    for( i = 0; i < p.size( ); i++ )
        p[ i ] *= rate;
}

#endif // OLD_STUFF

//=============================================================

void ChromosomeT< double >::initializeGSA( const ChromosomeT< double >&
					   sigma, int baseSize )
{

  initializeGSA( static_cast< const std::vector< double >& >( sigma ), 
			  baseSize );
  
}
//=============================================================
	     
void ChromosomeT< double >::initializeGSA( const std::vector< double >& 
					   sigma, int baseSize )
{

  int i,j,dim;
  dim = int( ((  *this ).size() - 1. ) / double( baseSize ) );
  
  (* this)[ 0 ] = sigma[ 0 ];

  for (j=0; j<baseSize; j++)
    for (i=0; i<dim; i++) 
      (* this)[ j*dim + i + 1 ] = sigma[ j*dim + i + 1];

}


//=============================================================
	     
void ChromosomeT< double >::initializeGSA( double SigmaMin, double SigmaMax,
					   int baseSize )
{

  int i,j,dim;
  dim = int( ((  *this ).size() - 1. ) / double( baseSize ) );
  
  (* this)[ 0 ] = Rng::uni( SigmaMin, SigmaMax);
  
  for (j=0; j<baseSize; j++)
    for (i=0; i<dim; i++) {
      if (j) (* this)[ j*dim + i + 1 ] = Rng::gauss(0, 1./dim );
      else (* this)[ j*dim + i + 1 ] = 0;
    }
}

//=============================================================

void ChromosomeT< double >::mutateGSA( ChromosomeT<double> & sigma )
{
  int Dimension   = ( *this ).size();
  int baseSize    = int( (sigma.size()-1)/(Dimension*1.0) );
  double beta     = 1. / sqrt( (double)Dimension );
  double cu       = sqrt( (2.-beta)/beta );
  double cm       = 1. / sqrt( (double)baseSize ) * ( 1. + 1./baseSize );
  double xi_const = 1.5;

  mutateGSA( sigma, beta, xi_const, cu, cm );

}

//=============================================================

void ChromosomeT< double >::mutateGSA( ChromosomeT<double> & sigma,
						double beta,
						double xi_const, double cu,
						double cm )
{
  
  int Dimension   = ( *this ).size();
  int baseSize    = int( (sigma.size()-1)/(Dimension*1.0) );
  double xi, *zet, **be, *ueps;
  int i,j;
  
  zet  = new double  [ baseSize ];
  ueps = new double  [ Dimension ];
  be   = new double* [ baseSize ];
  for (i=0; i < baseSize; i++) be[ i ] = new double [ Dimension ];
  
  if ( Rng::coinToss( 0.5 ) ) xi = xi_const;
  else xi = 1. / xi_const;
  
  // set the random vector z and for convenience copy the
  // sigma values from dad to a proper matrix B
  for (j=0; j < baseSize; j++ ) {
    zet[ j ] = Rng::gauss(0,1);
    for (i=0; i < Dimension; i++ )
      be[ j ][ i ] = sigma[ j * Dimension + i + 1 ];
  }
  
  for (i=0; i<Dimension; i++) {
    for (j=0, ueps[ i ] = 0; j < baseSize; j++ )
      ueps[ i ] += cm * ( zet[ j ] * be[ j ][ i ] );
  }

  // adaptation of the object variable
  for (i = 0; i<Dimension; i++) {
    ( *this )[ i ] += sigma[ 0 ] * xi * ueps[ i ];
  }
  
  // adaptation of the strategy variable
  // the global parameter firscht
  sigma[ 0 ]    = sigma[ 0 ] * pow( xi, beta );
  // the generating set b1,...,bm
  for (j=0; j < baseSize; j++ ) {
    for (i=0; i<Dimension; i++)
      if ( j ) sigma[ j * Dimension + i + 1 ] = be[ j-1 ][ i ];
      else sigma[ i + 1 ] = (1 - beta)*be[ 0 ][ i ] +
	     beta*cu*xi*ueps[ i ];
  }

  delete[ ] zet;
  delete[ ] ueps;
  for (i=0; i < baseSize; i++) delete [ ] be[ i ]; 
  delete [ ] be;

}  

//=============================================================

void ChromosomeT< double >::initializeIDA( const ChromosomeT< double >&
						    sigma )
{

  initializeIDA( static_cast< const std::vector< double >& >( sigma ) );
  
}
//=============================================================
	     
void ChromosomeT< double >::initializeIDA( const std::vector< double >& 
						    sigma )
{

  int i;
  
  for (i=0; i< int(sigma.size()); i++) (* this)[ i ] = sigma[ i ];

}

//=============================================================

void ChromosomeT< double >::initializeIDA( double SigmaMin, 
						    double SigmaMax )
{
  
  int i, dim;

  dim = int( ( (* this).size() - 2)/3. );
  // step size delta_r
  (* this)[ 0 ] = Rng::uni( SigmaMin, SigmaMax );
  // cumulation direction
  (* this)[ 2*dim + 1 ] = 0;
  // individual step size delta_i 
  for (i = 0; i<dim; i++) (* this)[ i + 1 ] = Rng::uni( SigmaMin, SigmaMax );
  // direction vector
  for (i = 0; i<dim; i++) (* this)[ dim + i + 1 ];
  // accumulated individual steps
  for (i = 0; i<dim; i++) (* this)[ 2*dim + i + 2 ];
}

//=============================================================

//#ifndef _WIN32 // fatal error C1001: INTERNER COMPILER-FEHLER

void ChromosomeT< double >::mutateIDA( ChromosomeT<double> & sigma ) {

  int Dimension   = ( *this ).size();
  double c        = 1. / sqrt( (double) Dimension );
  double c_r      = 3. / Dimension;
  double beta     = 2. / Dimension;
  double beta_ind = 1. / (4.*Dimension);
  double beta_r   = 1. / sqrt(4.*Dimension);
  double cu       = sqrt( (2.-c)/c );
  double xi       = 1/3.;

  mutateIDA( sigma, c, c_r, beta, beta_ind, beta_r, cu, xi );

}

//=============================================================

void ChromosomeT< double >::mutateIDA( ChromosomeT<double> & sigma, 
						double c, 
						double c_r, double beta, 
						double beta_ind, double beta_r, 
						double cu, double xi ) {
  
  int Dimension   = ( *this ).size();
  double chi_n    = sqrt( (double)Dimension )*( 1 - 1. / (4.*Dimension)
					      +  1. / (21.*sqr((double)Dimension)) );
  double chi_1    = sqrt( 2. / 3.141592654 );
    
  double aux, delta_r, s_r, zet_r, s_rN, delta_rN;
  double Norm_sN_vec, Norm_rN_vec, Norm_deltaN_vec;
  double *delta, *deltaN, *r, *rN, *s, *sN, *zet, *dist;

  int i;
  
  delta  = new double[ Dimension ];
  deltaN = new double[ Dimension ];
  s      = new double[ Dimension ];
  sN     = new double[ Dimension ];
  r      = new double[ Dimension ];
  rN     = new double[ Dimension ];
  zet    = new double[ Dimension ];
  dist   = new double[ Dimension ];

  // this stuff here is very close to (Hansen, 1995), the sections and
  // equations are given whenever appropriate, we do not give them if
  // it is not appropriate because that would be a waste of time and
  // we don't wanna do that now , do we...
  // extract the relevant variables from sigma
  delta_r = sigma[ 0 ];
  s_r     = sigma[ 2*Dimension + 1 ] ;
  for (i = 0; i<Dimension; i++) delta[ i ] = sigma[ i + 1 ];
  for (i = 0; i<Dimension; i++) r[ i ]     = sigma[ Dimension + i + 1 ];
  for (i = 0; i<Dimension; i++) s[ i ]     = sigma[ 2*Dimension + i + 2 ];

  // everything marked with an N stands for N as in offspriNg :-)
  // adaptation of the object variable (2.5.2 - 1.)
  zet_r   = Rng::gauss(0,1);
  for (i = 0; i<Dimension; i++) {
    zet   [ i ]  =  Rng::gauss(0,1);
    dist  [ i ]  = delta[ i ] * zet[ i ] + zet_r * delta_r * r[ i ];
    ( *this )[ i ] += dist[ i ];
  }

  // adaptation of the individual step sizes (2.5.2 - 2.1.)
  for (i = 0, Norm_sN_vec = 0; i<Dimension; i++) {
    sN[ i ]     = (1.-c) * s[ i ] + c * cu * zet[ i ];
    Norm_sN_vec += sqr(sN[ i ]);
  }
  Norm_sN_vec = sqrt( Norm_sN_vec );
  // adaptation of the individual step sizes (2.5.2 - 2.2.)
  for (i = 0, Norm_deltaN_vec = 0; i<Dimension; i++) {
    deltaN[ i ] = delta[ i ] * exp( beta * ( Norm_sN_vec - chi_n) ) *  
      exp( beta_ind * ( fabs( sN[ i ] ) - chi_1  ) );
    Norm_deltaN_vec += sqr( deltaN[ i ] );
  }
  Norm_deltaN_vec = sqrt( Norm_deltaN_vec );
  
  // direction adaptation (2.5.2 - 3.1.)
  aux = (1.-c) * s_r + c * cu * zet_r;
  if (aux > 0) s_rN = aux; else s_rN = 0;
  // direction adaptation (2.5.2 - 3.2.)
  for (i = 0, Norm_rN_vec = 0; i<Dimension; i++) {
    rN[ i ]     = (1.-c_r) * r[ i ] * delta_r
      + c_r * dist[ i ];
    Norm_rN_vec += sqr(rN[ i ]);
  }
  Norm_rN_vec = sqrt( Norm_rN_vec );
  // direction adaptation (2.5.2 - 3.3.)
  for (i = 0; i<Dimension; i++) rN[ i ]     = rN[ i ] / Norm_rN_vec;
  // direction adaptation (2.5.2 - 3.4.)
  aux = delta_r * exp( beta_r * ( fabs( s_rN ) - chi_1 ) );
  if (aux > (Norm_deltaN_vec * xi) ) delta_rN = aux;
  else delta_rN = Norm_deltaN_vec/3.;
  
  // and now we write all the little Ns back into the big, big
  // sigma chromosome ...
  sigma[ 0 ] = delta_rN;
  sigma[ 2*Dimension + 1 ] = s_rN ;
  for (i = 0; i<Dimension; i++) sigma[ i + 1 ] = deltaN[ i ];
  for (i = 0; i<Dimension; i++) sigma[ Dimension + i + 1 ] = rN[ i ];
  for (i = 0; i<Dimension; i++) sigma[ 2*Dimension + i + 2 ] = sN[ i ];

  delete[ ] zet;
  delete[ ] deltaN;
  delete[ ] delta;
  delete[ ] s;
  delete[ ] sN;
  delete[ ] r;
  delete[ ] rN;
  delete[ ] dist;
}
//#endif

//=============================================================

void ChromosomeT< double >::mutateRotate( ChromosomeT<double> & sigma ) 
{
  int    Dimension   = ( *this ).size();  
  double epsi       = 0.0001;
  int    sigmaCheck = 0;
  double tau1       = 1. / sqrt( 2.* Dimension );
  double tau2       = 1. / sqrt( 2.* sqrt( Dimension*1.0 ) );
  double beta       = 0.0873;

  mutateRotate( sigma, tau1, tau2, beta, sigmaCheck, epsi );
  
}
  
//=============================================================

void ChromosomeT< double >::mutateRotate( ChromosomeT<double> & sigma_sqr,
					  double tau1, double tau2,
					  double beta, int sigmaCheck,
					  double epsi)

{

  int    Dimension   = ( *this ).size();  
  double alpha = 0, *s, *ss, rdno;
  
  int i, j, k, index = 0, sign_of_sigma;

  s       = new double[ Dimension ];
  ss	  = new double[ Dimension ];

  // here the adaptation parameters are mutated first!
  rdno = Rng::gauss(0,1);
  for (i = 0; i<Dimension; i++) {
    sigma_sqr[ i ] *= exp( tau1 * rdno + tau2 * Rng::gauss(0,1) );
    // check that sigma is not too small!
    if ( (sigmaCheck) && (sigma_sqr[ i ] < epsi*fabs(( *this )[ i ])) )
      sigma_sqr[ i ] = epsi*fabs(( *this )[ i ]);
  }
  for (i = 0; i<Dimension-1; i++) {
    for (j = i+1; j<Dimension; j++) {
      index = int( Dimension*(i+1) -i-1+j -0.5*i*(i+1) );
      sigma_sqr[ index ] += beta * Rng::gauss(0,1);
      // check that angles are within limits!
      if ( fabs( sigma_sqr[ index ] ) > 3.141592654 ) {
	if (sigma_sqr[ index ]>=0) sign_of_sigma=1; else sign_of_sigma=-1;
	sigma_sqr[ index ] = sigma_sqr[ index ] - 2 * 3.141592654
	  * sign_of_sigma;
      }
    }
  }
  
  // draw rd. numbers according to sigma values!
  for (k = 0; k< Dimension; k++)  
  {
    s[ k ] = Rng::gauss(0, sigma_sqr[ k ] );
    ss[k]  = 0;
  }

  // do the funny rotation bits!
  for (i = 0; i<Dimension-1; i++) {
    for (j = i+1; j<Dimension; j++) {
      for (k = 0; k< Dimension; k++) {
	alpha = sigma_sqr[ int( Dimension*(i+1) -i-1+j -0.5*i*(i+1) ) ];
	if (i==k) 
	    ss[ k ] = s[ k ] * cos(alpha ) - s[ j ] * sin( alpha );
	else 
	if (j==k) 
	    ss[ k ] = s[ i ] * sin( alpha )+s[ k ] * cos( alpha );
        else 
	    ss[k] = s[k];
      }
      for(k=0; k<Dimension; k++) s[k]=ss[k];
    }
  }
  
  // change the objective parameters
  for (i = 0; i< Dimension; i++)
    ( *this )[ i ] += s[ i ];
  
  delete[ ] s;
  delete[ ] ss;
}

//=============================================================

void ChromosomeT< double >::initializeRotate( double SigmaMin,
					      double SigmaMax )      
{
  int i, j;
  int Dimension = int( -0.5 + sqrt( 0.25 + 2*( *this ).size()) );
  
  for (i = 0; i<Dimension; i++)
    ( *this )[ i ] = Rng::uni( SigmaMin, SigmaMax );;
  
  for (i = 0; i<Dimension-1; i++) 
    for (j = i+1; j<Dimension; j++) 
      ( *this )[ int( Dimension*(i+1) -i-1+j -0.5*i*(i+1) ) ] = 0;
      
}

//=============================================================

void ChromosomeT< double >::initializeRotate( const ChromosomeT< double >&
					      sigma )
{

  initializeIDA( static_cast< const std::vector< double >& >( sigma ) );
  
}
//=============================================================
	     
void ChromosomeT< double >::initializeRotate( const std::vector< double >& 
					      sigma )
{

  int i;
  
  for (i=0; i< int(sigma.size()); i++) (* this)[ i ] = sigma[ i ];

}

//=============================================================

void ChromosomeT< double >::mutateCMA( ChromosomeT<double> & sigma )
{
  int Dimension   = (* this).size();
  double c        = 1. / sqrt( (double)Dimension );
  double cu       = sqrt( (2.-c) * c );
  double beta     = 1. / Dimension;
  double ccov     = 2. / sqr( (double)Dimension );

  mutateCMA( sigma, c, cu, ccov, beta );
  
}

//=============================================================

void ChromosomeT< double >::mutateCMA( ChromosomeT<double> & sigma,
					     double c, double cu,
					     double ccov, double beta )
{
  
  int Dimension   = (* this).size();
  double chi_n    = sqrt( (double)Dimension )*( 1 - 1. / (4.*Dimension)
					+  1. / (21.*sqr((double)Dimension)) );
  
  double delta, *s, *zet, *sN, deltaN, norm, normS;
  double *sDelta, *sDeltaN;
  
  int i, j;

  Array< double > B         ( Dimension, Dimension );
  Array< double > Bdelta    ( Dimension, Dimension );
  Array< double > C         ( Dimension, Dimension );
  Array< double > lambda    ( Dimension );
  Array< double > prod      ( Dimension );
  Array< double > prodDelta ( Dimension );

  s       = new double[ Dimension ];
  zet     = new double[ Dimension ];
  sN      = new double[ Dimension ];
  sDelta  = new double[ Dimension ];
  sDeltaN = new double[ Dimension ];
  
  delta = sigma[ 0 ];
  for (i = 0; i<Dimension; i++) {
    zet   [ i ] = Rng::gauss(0,1);
    s     [ i ] = sigma[ i + 1];
    sDelta[ i ] = sigma[ Dimension + i + 1];
    for (j = 0; j<Dimension; j++)
      C(i, j) = sigma[ (i+2)*Dimension + j + 1 ];
  }

  eigensymm( C, B, lambda );
  
  for (i = 0; i<Dimension; i++) {
    for (j = 0, norm = 0; j<Dimension; j++) norm += sqr( B(j, i) );
    norm = sqrt( norm );
    for (j = 0; j<Dimension; j++) {
      Bdelta(j, i) = B(j, i) / norm;
      B     (j, i) = sqrt( fabs(lambda( i ))) * Bdelta(j, i);
    }
  }
    
  // change the objective parameters
  for (i = 0; i<Dimension; i++) {
    for (j = 0, prod( i ) = 0, prodDelta( i ) = 0; j<Dimension; j++) {
      prod     ( i ) += B     (i,j) * zet[ j ];
      prodDelta( i ) += Bdelta(i,j) * zet[ j ];
    }
    (* this)[ i ] += delta * prod( i );
  }

  // change the adaptation parameters
  for (i = 0, normS = 0; i<Dimension; i++) {
    sN[ i ]      = (1 - c) * s[ i ]      + cu * prod( i );
    sDeltaN[ i ] = (1 - c) * sDelta[ i ] + cu * prodDelta( i );
    normS += sqr( sDeltaN[ i ] );
  }
  normS  = sqrt( normS );
  deltaN = delta * exp( beta * ( normS - chi_n ) );
  // change the covariance matrix
  for (i = 0; i<Dimension; i++) 
    for (j = 0; j<Dimension; j++) 
      C(i, j ) = ( 1 - ccov ) * C(i, j) + ccov * sN[ i ] * sN[ j ];
  
  // write the stuff back!
  sigma[ 0 ] = deltaN;
  for (i = 0; i<Dimension; i++) {
    sigma[ i + 1 ] = sN[ i ];
    sigma[ i + 1 + Dimension] = sDeltaN[ i ];
    for (j = 0; j<Dimension; j++)
      sigma[ (i+2)*Dimension + j + 1 ] = C(i, j);
  }

  delete[ ] zet;
  delete[ ] sDelta;
  delete[ ] sDeltaN;
  delete[ ] s;
  delete[ ] sN;
  
}

//=============================================================

void ChromosomeT< double >::initializeCMA( double SigmaMin,
					   double SigmaMax  )      
{
  int i, j;

  int Dimension = int ( -1 + sqrt( (double) (( *this ).size() )) );
  
  ( *this )[ 0 ] = Rng::uni( SigmaMin, SigmaMax );
  for (i = 0; i<Dimension; i++) {
    ( *this )[ i + 1 ] = 0;
    ( *this )[ Dimension + i + 1 ] = 0;
    for (j = 0; j<Dimension; j++)
      if (i!=j) ( *this )[ (i+2)*Dimension + j + 1 ] = 0;
      else ( *this )[ (i+2)*Dimension + j + 1 ] = 1;
  }
}

//=============================================================

void ChromosomeT< double >::initializeCMA( const ChromosomeT< double >&
					   sigma )
{

  initializeIDA( static_cast< const std::vector< double >& >( sigma ) );
  
}
//=============================================================
	     
void ChromosomeT< double >::initializeCMA( const std::vector< double >& 
					   sigma )
{

  int i;
  
  for (i=0; i< int(sigma.size()); i++) (* this)[ i ] = sigma[ i ];

}

//=============================================================

void ChromosomeT< double >::mutateMSR( double xi_prob )
{
  
  int    Dimension   = ( *this ).size();
  double xi;
  
  // here the adaptation parameters are mutated!
  for (int i = 0; i<Dimension; i++) {
    if (Rng::coinToss( xi_prob )) xi = 1.5;
    else xi = 1/1.5;
    ( *this )[ i ] *= xi;
  }
  
}

//=============================================================

void ChromosomeT< double >::showRotate( )
{

  int i, j;
  
  int n   = int( ( *this ).size() );
  int dim = int( - 0.5 + sqrt( 2.0 * n + 0.25 ) );
  
  std::cout << "sigma       = ";
  for (i = 0; i<dim; i++)
	  std::cout << ( *this )[ i ] << "  ";
  std::cout << std::endl;
  std::cout << "alpha       = ";
  for (i = 0; i<dim-1; i++) {
    for (j = i+1; j<dim; j++) {
		std::cout << ( *this )[ int( dim*(i+1) -i-1+j -0.5*i*(i+1) ) ] << "  ";
    }
	std::cout << std::endl << "             ";
  }
  std::cout << std::endl;

}

//=============================================================

void ChromosomeT< double >::showCMA( ) 
{

  int i, j;
  int n   = int( ( *this ).size() );
  int dim = int( sqrt( 1.0 * n ) - 1 );
  
  std::cout << "delta   = " << ( *this )[ 0 ] << std::endl;
  std::cout << "s       = ";
  for (i = 0; i<dim; i++)
	  std::cout << ( *this )[ i + 1 ] << "  ";
  std::cout << std::endl;
  std::cout << "sDelta  = ";
  for (i = 0; i<dim; i++)
	  std::cout << ( *this )[ i + dim + 1 ] << "  ";
  std::cout << std::endl;
  std::cout << "C       = ";
  for (i = 0; i<dim; i++) {
    for (j = 0; j<dim; j++) {
		std::cout << ( *this )[ (i+2)*dim + j + 1 ] << "  ";
    }
	std::cout << std::endl << "          ";
  }
  std::cout << std::endl;
  
}
//=============================================================

void ChromosomeT< double >::showIDA( ) 
{

  int i;
  int n = int( ( *this ).size() );
  
  std::cout << "delta_r = " << ( *this )[ 0 ] << std::endl;
  std::cout << "s_r     = " << ( *this )[ 2*n + 1 ] << std::endl;
  std::cout << "delta   = ";
  for (i = 0; i<n; i++)
	  std::cout << ( *this )[ i + 1 ] << "  ";
  std::cout << std::endl;
  std::cout << "r       = ";
  for (i = 0; i<n; i++)
	  std::cout << ( *this )[ n + i + 1 ] << "  ";
  std::cout << std::endl;
  std::cout << "s       = ";
  for (i = 0; i<n; i++)
	  std::cout << ( *this )[ 2 * n + i + 2 ] << "  ";
  std::cout << std::endl;
}

//=============================================================

void ChromosomeT< double >::showGSA( int baseSize) 
{

	int i, j, dim;
	std::cout << "global step size = " << ( *this )[ 0 ] << std::endl;
	dim = int( ((  *this ).size() - 1. ) / double( baseSize ) );
  
	for (j = 0; j < baseSize; j++)
	{
		std::cout << "base vector " << j << std::endl;
		for (i = 0; i < dim; i++) 
			std::cout << "( " << i << ", " <<
		( *this )[ dim * j + i + 1 ] << " )  ";
		std::cout << std::endl;
	}
}

//=============================================================

void ChromosomeT< double >::initializeIDAiso( const ChromosomeT< double >& sigma )
{

  initializeIDAiso( static_cast< const std::vector< double >& >( sigma ) );
  
}
//=============================================================
	     
void ChromosomeT< double >::initializeIDAiso( const std::vector< double >& sigma )
{

  int i, dim;
  dim = size() / 2;
  
  for (i = 0; i < int(sigma.size()); i++) (* this)[ i ] = sigma[ i ];
}

//=============================================================

void ChromosomeT< double >::initializeIDAiso( double SigmaMin, 
					      double SigmaMax )
{
  
  int i, dim;

  dim = size() / 2;

  // individual step size delta_i 
  for (i = 0; i<dim; i++) (* this)[ i ] = Rng::uni( SigmaMin, SigmaMax );
  // direction vector
  for (i = 0; i<dim; i++) (* this)[ dim + i ] = 0.;
}

//=============================================================

void ChromosomeT< double >::mutateIDAiso( ChromosomeT<double> & sigma ) {

  int Dimension   = ( *this ).size();
  double c        = 1. / sqrt( (double)Dimension );
  double beta     = 2. / Dimension;
  double beta_ind = 1. / (4.*Dimension);
  double cu       = sqrt( (2.-c)/c );

  mutateIDAiso( sigma, c, beta, beta_ind, cu );

}

//=============================================================

void ChromosomeT< double >::mutateIDAiso( ChromosomeT<double> & sigma, 
					  double c, 
					  double beta, 
					  double beta_ind, 
					  double cu ) {
  
  int Dimension   = size();
  double chi_n    = sqrt( (double)Dimension )*( 1 - 1. / (4. * Dimension)
					+  1. / (21. * sqr((double)Dimension)) );
  double chi_1    = sqrt( 2. / 3.141592654 );
    
  double Norm_sN_vec;
  double *delta, *deltaN, *r, *rN, *s, *sN, *zet, *dist;

  int i;
  
  delta  = new double[ Dimension ];
  deltaN = new double[ Dimension ];
  s      = new double[ Dimension ];
  sN     = new double[ Dimension ];
  r      = new double[ Dimension ];
  rN     = new double[ Dimension ];
  zet    = new double[ Dimension ];
  dist   = new double[ Dimension ];

  // this stuff here is very close to (Hansen, 1995), the sections and
  // equations are given whenever appropriate, we do not give them if
  // it is not appropriate because that would be a waste of time and
  // we don't wanna do that now , do we...
  // extract the relevant variables from sigma
  for (i = 0; i<Dimension; i++) delta[ i ] = sigma[ i  ];
  for (i = 0; i<Dimension; i++) s[ i ]     = sigma[ Dimension + i ];


  // everything marked with an N stands for N as in offspriNg :-)
  // adaptation of the object variable (2.5.2 - 1.)
  for (i = 0; i < Dimension; i++) {
    zet   [ i ]  = Rng::gauss(0,1);
    dist  [ i ]  = delta[ i ] * zet[ i ];
    ( *this )[ i ] += dist[ i ];
  }

  // adaptation of the individual step sizes (2.5.2 - 2.1.)
  for (i = 0, Norm_sN_vec = 0; i<Dimension; i++) {
    sN[ i ]     = (1.-c) * s[ i ] + c * cu * zet[ i ];
    Norm_sN_vec += sqr(sN[ i ]);
  }
  Norm_sN_vec = sqrt( Norm_sN_vec );
  // adaptation of the individual step sizes (2.5.2 - 2.2.)
  for (i = 0; i<Dimension; i++) {
    deltaN[ i ] = delta[ i ] * exp( beta * ( Norm_sN_vec - chi_n) ) *  
      exp( beta_ind * ( fabs( sN[ i ] ) - chi_1  ) );
  }
  
  
  // and now we write all the little Ns back into the big, big
  // sigma chromosome ...
  for (i = 0; i<Dimension; i++) sigma[ i ] = deltaN[ i ];
  for (i = 0; i<Dimension; i++) sigma[ Dimension + i ] = sN[ i ];

  delete[ ] zet;
  delete[ ] deltaN;
  delete[ ] delta;
  delete[ ] s;
  delete[ ] sN;
  delete[ ] r;
  delete[ ] rN;
  delete[ ] dist;
}

//===========================================================================

/*
#ifndef __NO_GENERIC_IOSTREAM

void ChromosomeT< double >::writeTo( std::ostream& os ) const
{
    os << "ChromosomeT<" << typeid( double ).name( ) << ">(" << size() << ")" << endl;
    for( unsigned i = 0; i < size( ); i++ ) {
        if( i ) os << '\t';
	os << ( *this )[ i ];
    }
    os << std::endl;
}

void ChromosomeT< double >::readFrom( std::istream& is )
{
    string s;
  //is.getline( s );
    is >> s;
    is.get( );   // skip end of line

    if( is.good( ) &&
	s.substr( 0, 12 ) == "ChromosomeT<" &&
	s.find( '>' ) != string::npos &&
	s.substr( 12, s.find( '>' ) - 12 ) == typeid( double ).name( ) ) {

        resize( atoi( s.substr( s.find( '>' ) + 2 ).c_str( ) ) );
	for( unsigned i = 0; i < size( ); i++ )
	    is >> ( *this )[ i ];
    } else
      is.setf( ios::failbit );
}

#endif // !__NO_GENERIC_IOSTREAM
*/

// for unconstrained problems
void ChromosomeT< double >::SBX(ChromosomeT< double >& mate, 
				double nc, double p) {
  unsigned i, n = (*this).size();
  double beta, x, u = 0.;

  if(n != mate.size()) {
    std::cerr << "SBX is only defined for chromosomes of equal length" 
	      << std::endl;
    exit(EXIT_FAILURE);
  }

  for(i = 0; i < n; i++) {
    if(Rng::coinToss(p)) {
      do { u = Rng::uni(0, 1); } while(u == 1.);
      if(u <= .5) 
	beta = pow(2 * u, 1. / (nc + 1.)); 
      else
	beta = pow(1. / (2 - 2 * u), 1. / (nc + 1.)); 
      x = ( *this )[i];
      ( *this )[i] = .5 * ( (1 + beta) * x + (1 - beta) * mate[i]);
      mate[i]    = .5 * ( (1 - beta) * x + (1 + beta) * mate[i]);
    }
  }
}

// for constrained problems
void ChromosomeT< double >::SBX(ChromosomeT< double >& mate , 
				double lower, 
				double upper, 
				double nc, double p,
				double epsilon) {
  unsigned i, n = (*this).size();
  double beta, betaQ, alpha, expp, y1=0, y2=0, u = 0.;

  if(n != mate.size()) {
    std::cerr << "SBX is only defined for chromosomes of equal length" 
	      << std::endl;
    exit(EXIT_FAILURE);
  }

  for(i = 0; i < n; i++) {
    if(Rng::coinToss(p)) {
      if(mate[i] < (*this)[i])
	{
	  y1 = mate[i]; 
	  y2 = (*this)[i];
	}
      else
	{
	  y1 = (*this)[i]; 
	  y2 = mate[i];
	}
      if(fabs(y2 - y1) > epsilon) 
	{ //  -> from Deb's implementation, not contained in any paper: prevents division by zero
	  // Find beta value
	  if( ( y1 -  lower ) > ( upper - y2 ) ) 
	    {
	      beta = 1 + (2*( upper - y2)/( y2 - y1 ));
	    } else 
	    {
	      beta = 1 + (2*( y1 -  lower )/( y2 - y1 ));
	    }

	  expp = (nc + 1.);
	  beta = 1. / beta;
	  
	  // Find alpha
	  alpha = 2. - pow( beta , expp) ;

	  if (alpha < 0.0) 
	    {
	      std::cerr << "Error: " << alpha << " " << (*this)[i] << " " <<  mate[i] << std::endl;
	      exit(EXIT_FAILURE);
	    }

	  expp = 1./expp;

	  u = Rng::uni(0, 1);
	  //  -> from Deb's implementation, not contained in any paper 
	  // do { u = Rng::uni(0, 1); } while(u == 1.);

	  if(u <= 1./alpha) 
	    {
	      alpha *= u;
	      betaQ = pow(alpha, expp); 
	    }
	  else
	    {
      	      alpha *= u;
	      alpha = 1./( 2. - alpha );
	      if (alpha < 0.0) 
		{
		  std::cerr << "Error: " << alpha << " " << (*this)[i] << " " <<  mate[i] << std::endl;
		  exit(EXIT_FAILURE);
		}
	      betaQ = pow(alpha, expp); 
	    }
	} 
      else
	{ // if genes are equal -> from Deb's implementation, not contained in any paper
	  betaQ = 1.;
	}
      
      ( *this )[i]   = .5 * ( (y1 + y2) - betaQ * (y2 - y1) );
      mate[i]        = .5 * ( (y1 + y2) + betaQ * (y2 - y1) );
 
      //  -> from Deb's implementation, not contained in any paper
      if ( ( *this )[i] < lower) ( *this )[i] = lower;
      if ( ( *this )[i] > upper) ( *this )[i] = upper;
      if ( mate[i]      < lower) mate[i]      = lower;
      if ( mate[i]      > upper) mate[i]      = upper;
    }
  }
}
// for constrained problems
void ChromosomeT< double >::SBX(ChromosomeT< double >& mate , 
				std::vector<double > & lower, 
				std::vector<double > & upper, 
				double nc, double p,
				double epsilon) {
  unsigned i, n = (*this).size();
  double beta, betaQ, alpha, expp, y1=0, y2=0, u = 0.;

  if(n != mate.size()) {
    std::cerr << "SBX is only defined for chromosomes of equal length" 
	      << std::endl;
    exit(EXIT_FAILURE);
  }

  for(i = 0; i < n; i++) {
    if(Rng::coinToss(p)) {
      if(mate[i] < (*this)[i])
	{
	  y1 = mate[i]; 
	  y2 = (*this)[i];
	}
      else
	{
	  y1 = (*this)[i]; 
	  y2 = mate[i];
	}
      if(fabs(y2 - y1) > epsilon) 
	{ //  -> from Deb's implementation, not contained in any paper: prevents division by zero

	  // Find beta value
	  if((y1 - lower[i]) > ( upper[i] - y2)) 
	    beta = 1 + (2*( upper[i] - y2)/(y2 - y1));
    	  else 
	    beta = 1 + (2*( y1 - lower[i] )/(y2 - y1));
	    

	  expp = (nc + 1.);
	  beta = 1. / beta;
	  
	  // Find alpha
	  alpha = 2. - pow( beta , expp) ;

	  if (alpha < 0.0) 
	    {
	      std::cerr << "Error: " << alpha << " " << (*this)[i] << " " <<  mate[i] << std::endl;
	      exit(EXIT_FAILURE);
	    }

	  expp = 1./expp;

	  u = Rng::uni(0, 1);

	  //  -> from Deb's implementation, not contained in any paper 
	  // do { u = Rng::uni(0, 1); } while(u == 1.);

	  if(u <= 1./alpha) 
	    {
	      alpha *= u;
	      betaQ = pow(alpha, expp); 
	    }
	  else
	    {
      	      alpha *= u;
	      alpha = 1./( 2. - alpha );
	      if (alpha < 0.0) 
		{
		  std::cerr << "Error: " << alpha << " " << (*this)[i] << " " <<  mate[i] << std::endl;
		  exit(EXIT_FAILURE);
		}
	      betaQ = pow(alpha, expp); 
	    }
	} 
      else
	{ // if genes are equal -> from Deb's implementation, not contained in any paper
	  betaQ = 1.;
	}
      
      ( *this )[i]   = .5 * ( (y1 + y2) - betaQ * (y2 - y1) );
      mate[i]        = .5 * ( (y1 + y2) + betaQ * (y2 - y1) );
 
      //  -> from Deb's implementation, not contained in any paper
      if ( ( *this )[i] < lower[i]) ( *this )[i] = lower[i];
      if ( ( *this )[i] > upper[i]) ( *this )[i] = upper[i];
      if ( mate[i]      < lower[i]) mate[i]      = lower[i];
      if ( mate[i]      > upper[i]) mate[i]      = upper[i];
    }
  }
}

// for unconstrained problems
void ChromosomeT< double >::simpleMutatePolynomial(double lower, double upper, 
						   double nm, double p) {
  unsigned i, n = (*this).size();
  double delta,  r = 0.;

  for(i = 0; i < n; i++) {
    if(Rng::coinToss(p)) {
      r = Rng::uni(0, 1);
      if(r < .5) 
	delta = pow(2. * r, 1. / (nm + 1.)) - 1;
      else
	delta = 1 - pow(2. - 2. * r, 1. / (nm + 1.));
      (*this)[i] += delta * (upper - lower);
    }
  }
}
// for unconstrained problems
void ChromosomeT< double >::simpleMutatePolynomial(std::vector<double > &lower, 
						   std::vector<double > &upper, 
						   double nm, double p) {
  unsigned i, n = (*this).size();
  double delta,  r = 0.;

  for(i = 0; i < n; i++) {
    if(Rng::coinToss(p)) {
      r = Rng::uni(0, 1);
      if(r < .5) 
	delta = pow(2 * r, 1. / (nm + 1.)) - 1;
      else
	delta = 1 - pow(2 - 2 * r, 1. / (nm + 1.));
      (*this)[i] += delta * (upper[i] - lower[i]);
    }
  }
}
// for constrained problems
void ChromosomeT< double >::mutatePolynomial(double lower, 
					     double upper, 
					     double nm, double p) 
{
  unsigned i, n = (*this).size();
  double delta, deltaQ, expp,  u = 0.;

  for(i = 0; i < n; i++) 
    {
      if(Rng::coinToss(p)) 
	{
	  u  = Rng::uni(0, 1);
       
	  if ( (*this)[i] <=  lower || (*this)[i] >= upper )
	    { //  -> from Deb's implementation, not contained in any paper
	      (*this)[i] = u * ( upper -  lower ) +  lower;
#ifdef DEBUG
	      std::cerr << "Warning: parameter out of bounds, random resetting ..." << std::endl;
#endif
	    }
	  else 
	    {
	      // Calculate delta
	      if( ( (*this)[i] -  lower ) < ( upper - (*this)[i] ) )
		delta = ( (*this)[i] -  lower ) / ( upper -  lower );
	      else
		delta = ( upper - (*this)[i] ) / ( upper -  lower );

	      delta = 1. - delta;
	      expp  = ( nm + 1. );
	      delta = pow( delta , expp );
	      expp  = 1. / expp;

	      if(u <= .5) 
		{	
		  deltaQ =  2. * u + ( 1 - 2.*u ) * delta;
		  deltaQ = pow( deltaQ, expp ) - 1. ;
		}
	      else
		{
		  deltaQ = 2. - 2. * u + 2. * ( u  - .5 ) * delta;
		  deltaQ = 1. - pow( deltaQ , expp );
		}
	      
	      (*this)[i] += deltaQ * ( upper -  lower );

	      //  -> from Deb's implementation, not contained in any paper
	      if ((*this)[i] < lower) (*this)[i] = lower; 
	      if ((*this)[i] > upper) (*this)[i] = upper;
	    }
	}
    }
}
// for constrained problems
void ChromosomeT< double >::mutatePolynomial(std::vector<double > &lower, 
					     std::vector<double > &upper, 
					     double nm, double p) 
{
  unsigned i, n = (*this).size();
  double delta, deltaQ, expp,  u = 0.;
 

  for(i = 0; i < n; i++) 
    {
  
      if(Rng::coinToss(p)) 
	{
	  u  = Rng::uni(0, 1);
	  if ( (*this)[i] <= lower[i] || (*this)[i] >= upper[i] )
	    { //  -> from Deb's implementation, not contained in any paper
	      (*this)[i] = u * ( upper[i] - lower[i] ) + lower[i];
#ifdef DEBUG
	      std::cerr << "Warning: parameter out of bounds, random resetting ..." << std::endl;
#endif
	    }
	  else 
	    {
	      // Calculate delta
	      if( ( (*this)[i] - lower[i] ) < ( upper[i] - (*this)[i] ) )
		delta = ( (*this)[i] - lower[i] ) / ( upper[i] - lower[i] );
	      else
		delta = ( upper[i] - (*this)[i] ) / ( upper[i] - lower[i] );

	      delta = 1. - delta;
	      expp  = ( nm + 1. );
	      delta = pow( delta , expp );
	      expp  = 1. / expp;

	      if(u <= .5) 
		{	
		  deltaQ =  2. * u + ( 1 - 2.*u ) * delta;
		  deltaQ = pow( deltaQ, expp ) - 1. ;
		}
	      else
		{
		  deltaQ = 2. - 2. * u + 2. * (u  - .5) * delta;
		  deltaQ = 1. - pow( deltaQ , expp );
		}
	      
	      (*this)[i] += deltaQ * (upper[i] - lower[i]);

	      //  -> from Deb's implementation, not contained in any paper
	      if ((*this)[i] < lower[i]) 
		(*this)[i]=lower[i]; 
		
	      if ((*this)[i] > upper[i])  
		(*this)[i]=upper[i];
	
	    }
	}
    }
}


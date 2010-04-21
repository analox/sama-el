/* ======================================================================
 *
 *  File    :  MOOTestFunctions.h
 *  Created :  2001-12-16
 *
 *  Description : Class TestFunction in MOO-EALib
 *  Author  : Tatsuya Okabe <tatsuya.okabe@honda-ri.de>
 *  Copyright (c) 2001-2004 GPL2 Honda Research Institute Europe GmbH
 * 
 *  \par Maintained by:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>

 *  \par Project:
 *      MOO-EALib
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: MOOTestFunctions.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: MOOTestFunctions.h,v $
 *      Revision 1.2  2006/04/30 17:14:09  christian_igel
 *      *** empty log message ***
 *
 *      Revision 1.1  2006/04/30 17:11:50  christian_igel
 *      more test functions
 *
 *      Revision 2.16  2005/01/31 12:30:59  shark-admin
 *      definition of M_PI was lost
 *
 *      Revision 2.15  2005/01/31 11:12:15  shark-admin
 *      std:: added
 *
 *      Revision 2.14  2005/01/31 09:15:33  shark-admin
 *      bug in approxEquidistantFrontDoIt() corrected
 *
 *      Revision 2.13  2005/01/25 17:28:39  shark-admin
 *      Warnings added, front sampling replaced by approximate equidistant front sampling
 *
 *      Revision 2.12  2005/01/14 15:00:52  shark-admin
 *      namespaces std:: added
 *
 *      Revision 2.11  2005/01/14 11:41:21  shark-admin
 *      methods for sampling objective values and some more for pareto fronts added
 *
 *      Revision 2.10  2005/01/05 16:53:37  shark-admin
 *      rotated problems updated
 *
 *      Revision 2.9  2005/01/04 16:40:56  shark-admin
 *      Sample Front added for DebRotated
 *
 *      Revision 2.6  2004/12/30 15:39:28  shark-admin
 *      constant M_PI=3.14159265358979323846 added for Win32
 *
 *      Revision 2.5  2004/12/30 14:14:44  shark-admin
 *      namespace info and #include <EALib/sqr.h> added
 *
 *      Revision 2.4  2004/12/30 13:25:32  shark-admin
 *      ZDT-Test functions and rotated MOO paraboloid added
 *
 *      Revision 2.3  2004/12/15 14:43:59  shark-admin
 *      compatibility with win32 ensured
 *
 *      Revision 2.2  2004/03/03 13:23:49  saviapbe
 *      Copyright information was changed.
 *
 *      Revision 2.1  2004/02/12 11:27:09  saviapbe
 *      Author's contact address, INI Header are modified.
 *
 *      Revision 2.0  2003/11/28 16:23:11  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *<BR>
 
 * ----------------------------------------------------------------------
 *
 *  This file is part of the MOO-EALib. This library is free software;
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
 * ====================================================================== */
//
// 	Authors message
//======================================================================
/*	Thank you very much for your interest to MOO-EALib.

	Since our company's name was changed on 1st, January, 2003,
	my E-mail address in the source codes were also changed.
	The current E-mail address (6th,Feb.,2004) is as follows:

	tatsuya.okabe@honda-ri.de.

	If you cannot contact me with the above E-mail address,
	you can also use the following E-mail address:

	t_okabe_de@hotmail.com.

	If you have any questions, please don't hesitate to
	ask me. It's my pleasure.

	Best Regards,
	Tatsuya Okabe

	*********************************************************
	Tatsuya Okabe
	Honda Research Institute Europe GmbH
	Carl-Legien-Strasse 30, 63073 Offenbach/Main, Germany
	Tel: +49-69-89011-745
	Fax: +49-69-89011-749
	**********************************************************/

////////////////////////////////////////////////////////////////////// MOO


#ifndef __MOOTESTFUNCTIONS_H
#define __MOOTESTFUNCTIONS_H

#define _USE_MATH_DEFINES
#include "ZNconfig.h"
#include <math.h>
#include "EALib/sqr.h"

#ifdef _WIN32
  #include <limits>
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif
#else
  #include <values.h>
#endif

#include <Array/ArraySort.h>
#include <Array/ArrayIo.h>


//
//************************************************************************
// Tools for equidistant test function sampling in objective space 
//************************************************************************
//

unsigned approxEquidistantFrontDoIt(Array<double> &v, Array<double> &w, unsigned n, double d) {
  double dCurrent, dPrevious; 
  unsigned i, k;

  k = 1;
  for(i = 1; i < v.dim(0) - 1; i++) {
    dCurrent =  sqrt(sqr(w(k - 1, 0) - v(i    , 0)) + sqr(w(k - 1, 1) - v(i    , 1)));
    dPrevious = sqrt(sqr(w(k - 1, 0) - v(i - 1, 0)) + sqr(w(k - 1, 1) - v(i - 1, 1)));
    if(dCurrent > d) { 
      if(dPrevious == 0) { // if the adjacent point is far, take it
	w[k] = v[i];
	k++;
      } else { // take current
	if(fabs(dCurrent - d) < fabs(dPrevious - d)) {
	  if(k >= n) { 
	    return k + 1;
	  }
	  w[k] = v[i];
	  k++;
	} else { // take previous
	  if(k >= n) { 
	    return k + 1;
	  }
	  w[k] = v[i - 1];
	  k++;
	  i--;
	}
      }
    } 
  }
  return k;
}

void approxEquidistantFront(Array<double> &v, Array<double> &w, unsigned n = 500, double eps = 0.000001) {
  double sum = 0;
  unsigned i, k;
  if(v.ndim() != 2) {
    std::cerr << "makeUpsilonFront: input Array has wrong dimension" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(v.dim(0) < n) {
    std::cerr << "makeUpsilonFront: input Array has not enough sample points" << std::endl;
    exit(EXIT_FAILURE);
  } 

  w.resize(n, 2, false);

  // sort w.r.t. first component
  sort2DBy1st(v);

  // compute length of front
  for(i =0; i < v.dim(0) - 1; i++) 
    if(v(i + 1, 0) - v(i, 0)) sum += (v(i + 1, 0) - v(i, 0)) * 
      sqrt(1 + sqr((v(i + 1, 1) - v(i, 1))  / (v(i + 1, 0) - v(i, 0))));

  // fix end
  w[0] = v[0];

  double dOld, d;
  dOld = d = sum / (n - 1.);

  k = approxEquidistantFrontDoIt(v, w,  n, d);
  while (k != n) {
    dOld = d;
    if(k < n) d *= (1. - eps);     
    if(k > n) d *= (1. + eps);     
    k = approxEquidistantFrontDoIt(v, w,  n, d);
    //cerr << k << " " << d << endl;
  } 

 std::cout << "d = " << d << std::endl;
 std::cerr << "d = " << d << std::endl;

  w[n - 1] = v[v.dim(0) - 1];

}

//
//************************************************************************
// Tools for equidistant test function sampling in objective space 
//************************************************************************
//


//
//************************************************************************
// Test function implementation of Tatsuya Okabe
//************************************************************************
//

//************************************************************************
// Sphere Test Function (Deb's SCH)
//************************************************************************
double SphereF1( const std::vector< double >& x )
{
	double    sum = 0.0;
	unsigned  n;

	n = x.size( );
	for ( unsigned i = n; i--; )
	{
		sum += x[ i ] * x[ i ];
	}

	return sum / (double)n;
}
//
double SphereF2( const std::vector< double >& x )
{
	double   sum = 0.0;
	unsigned n;

	n = x.size( );
	for ( unsigned i = n; i--; )
	{
		sum += ( x[ i ] - 2.0 ) * ( x[ i ] - 2.0 );
	}

	return sum / (double)n;
}

//************************************************************************
// Deb's Convex Test Function (ZDT 1)
//************************************************************************
double DebConvexF1( const std::vector< double >& x )
{
	unsigned n;
	n = x.size( );
	// Penalty
	for ( unsigned i = 0; i < n; i++ )
	{
		if ( x[ i ] > 1.0 || x[ i ] < 0.0 )
		{
			return 5.0;
		}
	}
	return x[ 0 ];
}
//
double DebConvexF2( const std::vector< double >& x )
{
	unsigned i;
	double   sum = 0.0;
	unsigned n;
	double   g, f1, f2;
	f1 = x[ 0 ];
	n  = x.size( );
	// Penalty
	for ( i = 0; i < n; i++ )
	{
		if ( x[ i ] > 1.0 || x[ i ] < 0.0 )
		{
			return 10.0;
		}
	}
	//
	for ( i = 1; i < n; i++ )
	{
		sum += x[ i ];
	}
	g = 1.0 + 9.0 * sum / (double)( n - 1 );
	f2 = g * ( 1.0 - sqrt( f1 / g ) );
	return f2;
}

//************************************************************************
// Deb's Concave Test Function (ZDT 2)
//************************************************************************
double DebConcaveF1( const std::vector< double >& x )
{
	unsigned n;
	unsigned i;
	n = x.size( );
	// Penalty
	for ( i = 0; i < n; i++ )
	{
		if ( x[ i ] > 1.0 || x[ i ] < 0.0 )
		{
			return 5.0;
		}
	}

	return x[ 0 ];
}
//
double DebConcaveF2( const std::vector< double >& x )
{
	double   sum = 0.0;
	unsigned n;
	double   g, f1, f2;
	unsigned i;
	f1 = x[ 0 ];
	n  = x.size( );
	// Penalty
	for ( i = 0; i < n; i++ )
	{
		if ( x[ i ] > 1.0 || x[ i ] < 0.0 )
		{
			return 10.0;
		}
	}
	//
	for ( i = 1; i < n; i++ )
	{
		sum += x[ i ];
	}
	g = 1.0 + 9.0 * sum / (double)( n - 1 );
	f2 = g * ( 1.0 - ( f1 / g ) * ( f1 / g ) );
	return f2;
}

//************************************************************************
// Deb's Discrete Test Function (ZDT 3)
//************************************************************************
double DebDiscreteF1( const std::vector< double >& x )
{
	unsigned n;
	unsigned i;
	n = x.size( );
	// Penalty
	for ( i = 0; i < n; i++ )
	{
		if ( x[ i ] > 1.0 || x[ i ] < 0.0 )
		{
			return 5.0;
		}
	}

	return x[ 0 ];
}
//
double DebDiscreteF2( const std::vector< double >& x )
{
	double   sum = 0.0;
	unsigned n;
	double   g, f1, f2;
	unsigned i;
	f1 = x[ 0 ];
	n  = x.size( );
	// Penalty
	for ( i = 0; i < n; i++ )
	{
		if ( x[ i ] > 1.0 || x[ i ] < 0.0 )
		{
			return 10.0;
		}
	}
	//
	for ( i = 1; i < n; i++ )
	{
		sum += x[ i ];
	}
	g = 1.0 + 9.0 * sum / (double)( n - 1 );
	f2 = g * ( 1.0 - sqrt( f1 / g ) - ( f1 / g ) * sin( 10 * znPiC * f1 ) );
	return f2;
}

//************************************************************************
// Fonseca's Concave Test Function (Deb's FON)
//************************************************************************
double FonsecaConcaveF1( const std::vector< double >& x )
{
	double   sum = 0.0;
	unsigned n;
	unsigned i;
	n = x.size( );
	for ( i = n; i--; )
	{
		sum += ( x[ i ] - 1.0 / sqrt( (double)n )) * ( x[ i ] - 1.0 / sqrt( (double)n ));
	}



	return 1.0 - exp( -sum );
}
//
double FonsecaConcaveF2( const std::vector< double >& x )
{
	double   sum = 0.0;
	unsigned n;
	unsigned i;
	n = x.size( );
	for ( i = n; i--; )
	{
		sum += ( x[ i ] + 1.0 / sqrt( (double)n )) * ( x[ i ] + 1.0 / sqrt( (double)n ));
	}

	return 1.0 - exp( -sum );
}

void FonsecaConcaveSampleFront(Array< double > &pf, unsigned dimension, unsigned n  = 500) 
{
        double xmin = -1./sqrt((double)dimension);
	double xmax = 1./sqrt((double)dimension);
	unsigned i  = 0, ii = 0;
	std::vector< double > xv(dimension);

	pf.resize(n, 2u);
	unsigned raw = 10 * n;
	Array<double> raw_pf(raw, 2u);

	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin;  x < xmax - 0.1 * it; x += it,i++) 
	{ 
		for(ii = 0; ii < xv.size(); ii++)
		  xv[ii] = x;
		
		raw_pf(i, 0) = FonsecaConcaveF1( xv );
		raw_pf(i, 1) = FonsecaConcaveF2( xv );
	}
	i = raw - 1;
	for(ii = 0; ii < xv.size(); ii++)
	{
	  xv[ii] = xmax;
	}
	raw_pf(i, 0) = FonsecaConcaveF1( xv );
	raw_pf(i, 1) = FonsecaConcaveF2( xv );

	approxEquidistantFront(raw_pf, pf, n);
}

void FonsecaConcaveSample(Array< double > &pf, unsigned dimension, unsigned n  = 500) 
{
	if(dimension!=2)
	{
		std::cerr << "sorry, method implemented for two dimensions only ..." << std::endl;
		return;
	}
	double xmin = -1./sqrt((double)dimension);
	double xmax = 1./sqrt((double)dimension);
	unsigned i  = 0;
	std::vector< double > xv(dimension);

	pf.resize(n*n, 2u);

	double it = (xmax - xmin) / (n - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it) 
	{
		xv[0] = x;
		for(double y = xmin; y < xmax - 0.1 * it; y += it) 
		{
			xv[1] = y; 
			pf(i, 0) = FonsecaConcaveF1( xv );
			pf(i, 1) = FonsecaConcaveF2( xv );
			i++;
		}
		xv[1] = xmax; 
		pf(i, 0) = FonsecaConcaveF1( xv );
		pf(i, 1) = FonsecaConcaveF2( xv );
		i++;
	}
	xv[0] = xmax;
	for(double y = xmin; y < xmax - 0.1 * it; y += it) 
	{
		xv[1] = y; 
		pf(i, 0) = FonsecaConcaveF1( xv );
		pf(i, 1) = FonsecaConcaveF2( xv );
		i++;
	}
	xv[1] = xmax; 
	pf(i, 0) = FonsecaConcaveF1( xv );
	pf(i, 1) = FonsecaConcaveF2( xv );
}

//************************************************************************
// Messac's Concave Test Function
//************************************************************************
double MessacConcaveF1( const std::vector< double >& x )
{
	double   sum = 0.0;
	unsigned n,i;

	n = x.size( );
	for ( i = n; i--; )
	{
		sum += exp( -x[ i ] ) + 1.4 * exp( -x[ i ] * x[ i ] );
	}

	return sum;
}
//
double MessacConcaveF2( const std::vector< double >& x )
{
	double   sum = 0.0;
	unsigned n,i;

	n = x.size( );
	for ( i = n; i--; )
	{
		sum += exp( +x[ i ] ) + 1.4 * exp( -x[ i ] * x[ i ] );
	}

	return sum;
}

void MessacConcaveSampleFront(Array< double > &pf, unsigned dimension, double lower, double upper, unsigned n  = 500) 
{
	if(dimension!=1)
	{
		std::cerr << "sorry, method implemented for one dimension only ..." << std::endl;
		return;
	}

	double xmin = lower;
	double xmax = upper;
	unsigned i = 0, ii=0;
	std::vector< double > xv(dimension);

	pf.resize(n, 2u);

	double it = (xmax - xmin) / (n - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{ 
		for(ii = 0; ii < xv.size(); ii++) xv[ii] = x;
		pf(i, 0) = MessacConcaveF1( xv );
		pf(i, 1) = MessacConcaveF2( xv );
	}
	i = n - 1;
	for(ii = 0; ii < xv.size(); ii++) xv[ii] = xmax;
	pf(i, 0) = MessacConcaveF1( xv );
	pf(i, 1) = MessacConcaveF2( xv );
}

void MessacConcaveSample(Array< double > &pf, unsigned dimension, double lower, double upper) 
{
	if(dimension!=2)
	{
		std::cerr << "sorry, method implemented for two dimensions only ..." << std::endl;
		return;
	}

	double xmin = lower;
	double xmax = upper;
	unsigned n  = 100, i = 0;
	std::vector< double > xv(dimension);

	pf.resize(n*n, 2u);

	double it = (xmax - xmin) / (n - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it) 
	{
		xv[0] = x;
		for(double y = xmin; y < xmax - 0.1 * it; y += it) 
		{
			xv[1] =  y;
			pf(i, 0) = MessacConcaveF1( xv );
			pf(i, 1) = MessacConcaveF2( xv );
			i++;
		}
		xv[1] =  xmax;
		pf(i, 0) = MessacConcaveF1( xv );
		pf(i, 1) = MessacConcaveF2( xv );
		i++;
	}
	xv[0] = xmax;
	for(double y = xmin; y < xmax - 0.1 * it; y += it) 
	{
		xv[1] =  y;
		pf(i, 0) = MessacConcaveF1( xv );
		pf(i, 1) = MessacConcaveF2( xv );
		i++;
	}
	xv[1] =  xmax;
	pf(i, 0) = MessacConcaveF1( xv );
	pf(i, 1) = MessacConcaveF2( xv );
}

//
//************************************************************************
// End of test function implementation of Tatsuya Okabe
//************************************************************************
//

//
//************************************************************************
// Rotated test functions (by Stefan Roth)
//************************************************************************
//

//************************************************************************
// Helpers
//************************************************************************

double norm(const Array<double> &a) 
{
	double sum = 0;
	unsigned i ;
	for(i = 0; i < a.nelem(); i++) sum += a(i)*a(i);
	return sqrt(sum);
}

double scalarprod(const Array<double> &a, const Array<double> &b) 
{
	double sum = 0;
	unsigned i;
	for(i = 0; i < a.nelem(); i++)
		sum += a(i) * b(i);
	return sum;
}

double scalarprod(const Array<double> &a, const std::vector<double> &b) 
{
	double sum = 0;
	unsigned i;

	if(a.ndim()!=1 || a.nelem()!=b.size())
	{
		std::cerr << "check size of vector or Array" << std::endl;
		exit(0);
	}

	for(i = 0; i < a.nelem(); i++)
		sum += a(i) * b[i];
	return sum;
}


void generateBasis(unsigned d, Array<double> &B) 
{
  unsigned i, j, c;
  double normB;
  Array<double> H;
  B.resize(d, d);
  H.resize(d, d);
  for(i = 0; i < d; i++) 
    for(c = 0; c < d; c++)
      B(i, c) = Rng::gauss(0, 1);

  for(i = 0; i < d; i++) {
    for(j = 0; j < i; j++) {
      H=B;
      for(c = 0; c < d; c++) 
	B(i, c) -= scalarprod(H[i], H[j]) * H(j, c);
    }
    normB = norm(B[i]);
    for(j = 0; j < d; j++) B(i, j) = B(i, j) / normB;
  }
}

//************************************************************************
// Rotated paraboloid 
//************************************************************************

double RotParF1(const std::vector<double> &_v, Array<double> &coord, double cond)
{
	unsigned i; 
	double sum = 0.;

	Array<double> v(_v.size());

	for(i = 0; i < v.dim(0); i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += sqr(pow(cond, (double(i) / double(v.dim(0) - 1))) * v(i)) ;
	}

	return sum / (cond * cond * v.dim(0));
}

double RotParF2(const std::vector<double> &_v, Array<double> &coord, double cond1, double cond2=2)
{
	unsigned i; 
	double sum = 0.;

	Array<double> v(_v.size());

	for(i = 0; i < v.dim(0); i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += sqr(pow(cond1, (double(i) / double(v.dim(0) - 1))) * (v(i) - cond2)) ;
	}

	return sum / (cond1 * cond1 * v.dim(0));
}

double RotParF1(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond)
{
  unsigned i; 

  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotParF1(_v, coord, cond);
}

double RotParF2(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond1, double cond2=2)
{
  unsigned i; 
  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotParF2(_v, coord, cond1,cond2);
}

void RotParSampleFront(Array< double > &pf, unsigned dimension, double cond , double cond2 = 2, unsigned n  = 100000) 
{
	double xmin = 0;
	double xmax = cond2;
	unsigned  i = 0,ii=0;
	std::vector< double > xv(dimension);

	pf.resize(n, 2u);
	unsigned raw = 40 * n;
	Array<double> raw_pf(raw, 2u);

	Array<double> B(dimension,dimension);
	B = 0; for(i=0;i<B.dim(0);i++) B(i,i)=1;

	i=0;
	double it = (xmax - xmin) / ( raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{ 
		for(ii = 0; ii < xv.size(); ii++) xv[ii] = x;
		raw_pf(i, 0) = RotParF1( xv, B, cond );
		raw_pf(i, 1) = RotParF2( xv, B, cond, cond2 );
	}
	i = raw - 1;
	for(ii = 0; ii < xv.size(); ii++) xv[ii] = xmax;
	raw_pf(i, 0) = RotParF1( xv, B, cond );
	raw_pf(i, 1) = RotParF2( xv, B, cond, cond2 );

	approxEquidistantFront(raw_pf, pf, n);
}

void RotParSample(Array< double > &pf, unsigned dimension, double cond , double cond2=2 ) 
{
	if(dimension!=2)
	{
		std::cerr << "sorry, method implemented for one dimension only ..." << std::endl;
		return;
	}

	double xmin = 0, x, y;
	double xmax = cond2;
	unsigned n  = 100, i = 0;
	std::vector< double > xv(dimension);

	pf.resize(n*n, 2u);

	Array<double> B(dimension,dimension);
	B = 0; for(i=0;i<B.dim(0);i++) B(i,i)=1;

	i = 0;

	double it = (xmax - xmin) / (n - 1.);

	for(x = xmin; x < xmax - 0.1 * it; x += it) 
	{
		xv[0] = x;
		for(y = xmin; y < xmax - 0.1 * it; y += it) 
		{
			xv[1] = y;
			pf(i, 0) = RotParF1( xv, B, cond );
			pf(i, 1) = RotParF2( xv, B, cond, cond2 );
			i++;
		}
		xv[1] = xmax;
		pf(i, 0) = RotParF1( xv, B, cond );
		pf(i, 1) = RotParF2( xv, B, cond, cond2 );
		i++;
	}
	xv[0] = xmax;
	for(y = xmin; y < xmax - 0.1 * it; y += it) 
	{
		xv[1] = y;
		pf(i, 0) = RotParF1( xv, B, cond );
		pf(i, 1) = RotParF2( xv, B, cond, cond2 );
		i++;
	}
	xv[1] = y;
	pf(i, 0) = RotParF1( xv, B, cond );
	pf(i, 1) = RotParF2( xv, B, cond, cond2 );
}


//************************************************************************
// Rotated cigar
//************************************************************************

double RotCigarF1(const std::vector<double> &_v, Array<double> &coord, double cond)
{
	unsigned i; 
	double sum = 0.;

	Array<double> v(_v.size());

	v(0) = scalarprod(coord.col(0),_v);
	sum += v(0) * v(0) ;

	for(i = 1; i < v.dim(0); i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += sqr( cond * v(i) ) ;
	}

	return sum / (cond * cond * v.dim(0));
}

double RotCigarF2(const std::vector<double> &_v, Array<double> &coord, double cond1, double cond2=2)
{
	unsigned i; 
	double sum = 0.;

	Array<double> v(_v.size());

	v(0) = scalarprod(coord.col(0),_v);
	sum += (v(0) - cond2) * (v(0) - cond2) ;

	for(i = 1; i < v.dim(0); i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += sqr( cond1 * (v(i) - cond2) ) ;
	}

	return sum / (cond1 * cond1 * v.dim(0));
}

double RotCigarF1(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond)
{
  unsigned i; 

  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotCigarF1(_v, coord, cond);
}

double RotCigarF2(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond1, double cond2=2)
{
  unsigned i; 
  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotCigarF2(_v, coord, cond1,cond2);
}

void RotCigarSampleFront(Array< double > &pf, unsigned dimension, double cond , double cond2=2, unsigned n  = 100000) 
{
	double xmin = 0;
	double xmax = cond2;
	unsigned i = 0, ii=0;
	std::vector< double > xv(dimension);

	pf.resize(n, 2u);
	unsigned raw = 10 * n;
	Array<double> raw_pf(raw, 2u);

	Array<double> B(dimension,dimension);
	B = 0; for(i=0;i<B.dim(0);i++) B(i,i)=1;

	i = 0;

	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{ 
		for(ii = 0; ii < xv.size(); ii++) xv[ii] = x;
		raw_pf(i, 0) = RotCigarF1( xv, B, cond );
		raw_pf(i, 1) = RotCigarF2( xv, B, cond, cond2 );
	}
	i = raw - 1;
	for(ii = 0; ii < xv.size(); ii++) xv[ii] = xmax;
	raw_pf(i, 0) = RotCigarF1( xv, B, cond );
	raw_pf(i, 1) = RotCigarF2( xv, B, cond, cond2 );

	approxEquidistantFront(raw_pf, pf, n);
}

//************************************************************************
// Rotated cigtab
//************************************************************************

double RotCigTabF1(const std::vector<double> &_v, Array<double> &coord, double cond)
{
	unsigned i, n; 
	double sum = 0.;

	n = _v.size();
	Array<double> v;
	v.resize(n, false);

	v(0) = scalarprod(coord.col(0),_v);
	sum += v(0) * v(0) ;
	v(n - 1) = scalarprod(coord.col(n - 1),_v);
	sum += sqr(cond * v(n - 1)) ;

	for(i = 1; i < n - 1; i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += cond * sqr( v(i) ) ;
	}

	return sum / (cond * cond * n);
}

double RotCigTabF2(const std::vector<double> &_v, Array<double> &coord, double cond1, double cond2=2)
{
	unsigned i, n; 
	double sum = 0.;

	n = _v.size();
	Array<double> v;
	v.resize(n, false);

	v(0) = scalarprod(coord.col(0),_v);
	sum += sqr(v(0) - cond2);
	v(n - 1) = scalarprod(coord.col(n - 1),_v);
	sum += sqr(cond1 * (v(n - 1) - cond2));

	for(i = 1; i < n - 1; i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += cond1 * sqr(v(i) - cond2) ;
	}

	return sum / (cond1 * cond1 * n);
}

double RotCigTabF1(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond)
{
  unsigned i; 

  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotCigTabF1(_v, coord, cond);
}

double RotCigTabF2(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond1, double cond2=2)
{
  unsigned i; 
  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotCigTabF2(_v, coord, cond1,cond2);
}

void RotCigTabSampleFront(Array< double > &pf, unsigned dimension, double cond , double cond2=2, unsigned n  = 100000) 
{
	double xmin = 0;
	double xmax = cond2;
	unsigned i = 0, ii=0;
	std::vector< double > xv(dimension);

	pf.resize(n, 2u);
	unsigned raw = 10 * n;
	Array<double> raw_pf(raw, 2u);

	Array<double> B(dimension,dimension);
	B = 0; for(i=0;i<B.dim(0);i++) B(i,i)=1;

	i = 0;

	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{ 
		for(ii = 0; ii < xv.size(); ii++) xv[ii] = x;
		raw_pf(i, 0) = RotCigTabF1( xv, B, cond );
		raw_pf(i, 1) = RotCigTabF2( xv, B, cond, cond2 );
	}
	i = raw - 1;
	for(ii = 0; ii < xv.size(); ii++) xv[ii] = xmax;
	raw_pf(i, 0) = RotCigTabF1( xv, B, cond );
	raw_pf(i, 1) = RotCigTabF2( xv, B, cond, cond2 );

	approxEquidistantFront(raw_pf, pf, n);
}

//************************************************************************
// Rotated tablet
//************************************************************************

double RotTabletF1(const std::vector<double> &_v, Array<double> &coord, double cond)
{
	unsigned i; 
	double sum = 0.;

	Array<double> v(_v.size());

	v(0) = scalarprod(coord.col(0),_v);
	sum += sqr(cond * v(0));

	for(i = 1; i < v.dim(0); i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += v(i) * v(i);
	}

	return sum / (cond * cond * v.dim(0));
}

double RotTabletF2(const std::vector<double> &_v, Array<double> &coord, double cond1, double cond2=2)
{
	unsigned i; 
	double sum = 0.;

	Array<double> v(_v.size());

	v(0) = scalarprod(coord.col(0),_v);
	sum += sqr(cond1 * (v(0) - cond2)) ;

	for(i = 1; i < v.dim(0); i++)
	{
		v(i) = scalarprod(coord.col(i),_v);
		sum += (v(i) - cond2) * (v(i) - cond2) ;
	}

	return sum / (cond1 * cond1 * v.dim(0));
}

double RotTabletF1(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond)
{
  unsigned i; 

  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotTabletF1(_v, coord, cond);
}

double RotTabletF2(const std::vector<double> &_v, Array<double> &coord, const std::vector<double> &lower, const std::vector<double> &upper, double cond1, double cond2=2)
{
  unsigned i; 
  for (i = 0; i < _v.size();i++)
    {
      // Penalty
      if ( _v[i] > upper[i] || _v[i] < lower[i] )
	{
#ifdef _WIN32
	  return std::numeric_limits< double >::max( );
#else
	  return MAXDOUBLE;
#endif
	}
      //
    }

  return RotTabletF2(_v, coord, cond1,cond2);
}

void RotTabletSampleFront(Array< double > &pf, unsigned dimension, double cond , double cond2 = 2, unsigned  n  = 100000) 
{
	double xmin = 0;
	double xmax = cond2;
	unsigned i = 0,ii=0;
	std::vector< double > xv(dimension);

	pf.resize(n, 2u);
	unsigned raw = 20 * n;
	Array<double> raw_pf(raw, 2u);

	Array<double> B(dimension,dimension);
	B = 0; for(i=0;i<B.dim(0);i++) B(i,i)=1;

	i = 0;
	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{ 
		for(ii = 0; ii < dimension; ii++) xv[ii] = x;
		raw_pf(i, 0) = RotTabletF1( xv, B, cond );
		raw_pf(i, 1) = RotTabletF2( xv, B, cond, cond2 );
	}
	i = raw - 1;
	for(ii = 0; ii < dimension; ii++) xv[ii] = xmax;
	raw_pf (i, 0) = RotTabletF1( xv, B, cond );
	raw_pf (i, 1) = RotTabletF2( xv, B, cond, cond2 );

	approxEquidistantFront(raw_pf, pf, n);
}

//************************************************************************
// Fonseca's Concave Test Function (Deb's FON)
//************************************************************************
double RotFonsecaConcaveF1( const std::vector< double >& x, Array<double> &coord, double cond)
{
	unsigned n = x.size();
	unsigned i;
	Array<double> y(n);
	for(i = 0; i < n; i++) y(i) = pow(cond, double(i) / double(n - 1)) * scalarprod(coord.col(i), x);
	double   sum = 0.0;

	for ( i = n; i--; )
	  sum += ( y(i) - 1.0 / sqrt( (double)n )) * ( y(i) - 1.0 / sqrt( (double)n ));

	return 1.0 - exp( -sum );
}
//
double RotFonsecaConcaveF2( const std::vector< double >& x, Array<double> &coord, double cond)
{
	unsigned n = x.size();
	unsigned i;
	Array<double> y(n);
	for(i = 0; i < n; i++) y(i) = pow(cond, double(i) / double(n - 1)) * scalarprod(coord.col(i), x);
	double   sum = 0.0;

	for ( i = n; i--; ) sum += ( y(i) + 1.0 / sqrt( (double)n )) * ( y(i) + 1.0 / sqrt( (double)n ));

	return 1.0 - exp( -sum );
}


//************************************************************************
// Deb's Rotated problem (IEEE Transactions on EA 6(2), 2002)
//************************************************************************

double DebRotatedF1( const std::vector< double > &x, Array<double> &coord )
{
	double f1;//, penalty = 0;
	f1 = scalarprod(coord.col(0),x);

	// Penalty
	if ( f1 > .3 || f1 < -.3 )
	{
		//penalty = fabs(f1) - .3;
#ifdef _WIN32
		return std::numeric_limits< double >::max( );
#else
		return MAXDOUBLE;
#endif
	}
	//
	return f1 ;//+ penalty;
}


//
double DebRotatedF2( const std::vector< double >& x, Array<double> &coord )
{
	double   sum = 0.0;
	unsigned n;
	double   g, f1, f2;//, penalty=0;
	unsigned i;

	n  = x.size( );

	Array<double> v(n);

	for(i = 0; i < n; i++)
		v(i) = scalarprod(coord.col(i),x);

	f1 = v(0);

	if ( f1 > .3 || f1 < -.3 )
	{
		//penalty = fabs(f1) - .3;		
#ifdef _WIN32
		return std::numeric_limits< double >::max( );
#else
		return MAXDOUBLE;
#endif
	}

	//
	for ( i = 1; i < n; i++ )
		sum += v( i ) * v( i ) - 10. * cos(4 * M_PI * v(i));

	g = 1.0 + 10. * (double)( n - 1 ) + sum;

	f2 = g * exp( - 1. * ( f1 / g ) );

	return f2 ;//+ penalty;
}


void DebRotatedSampleFront(Array< double > &pf, unsigned dimension, unsigned n  = 500) 
{
	double xmin = -.3;
	double xmax = .3;
	unsigned i = 0;
	std::vector< double > xv(dimension);

	pf.resize(n, 2u);
	unsigned raw = 10 * n;
	Array<double> raw_pf(raw, 2u);

	Array<double> B(dimension,dimension);
	B = 0; for(i=0;i<B.dim(0);i++) B(i,i)=1;

	for(i = 1; i < xv.size(); i++) xv[i] = 0.;

	i = 0;
	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{
		xv[0] = x;
		raw_pf(i, 0) = DebRotatedF1( xv, B );
		raw_pf(i, 1) = DebRotatedF2( xv, B );
	}
	i = raw - 1;
	xv[0] = xmax;
	raw_pf(i, 0) = DebRotatedF1( xv, B );
	raw_pf(i, 1) = DebRotatedF2( xv, B );

	approxEquidistantFront(raw_pf, pf, n);
}

//
//************************************************************************
// End of paraboloid test functions 
//************************************************************************
//

//
//************************************************************************
// ZDT test function implementation of Christian Igel
//************************************************************************
//

double CI1F1( const std::vector< double >& x )
{
	return x[0] * x[0];
}

double CI1F2( const std::vector< double >& x )
{
	return (x[1] - 2.) * (x[1] - 2.) + (x[2] - 2.) * (x[2] - 2.);
}

double CI2F1( const std::vector< double >& x )
{
	return x[0] * x[0];
}

double CI2F2( const std::vector< double >& x )
{
	return (x[0] - 2.) * (x[0] - 2.) + x[1] * x[1];
}

//************************************************************************
// ZDT 1
//************************************************************************

double ZDT1F1( const std::vector< double >& x )
{
	return fabs(x[0]);
}

double ZDT1G(const std::vector< double >& x )
{
	double g = 0;
	unsigned n = x.size( );
	for ( unsigned i = 1; i < n; i++ ) 
		g += fabs((x[i]));
	g = 9 * g / (n - 1.) + 1.;
	return g;
}

double ZDT1F2( const std::vector< double >& x )
{
	return  ZDT1G(x) * (1 - sqrt( fabs(x[0]) /  ZDT1G(x) ));
}

void ZDT1SampleFront(Array< double > &pf, unsigned n = 500) 
{
	double xmin = 0;
	double xmax = 1.;

	pf.resize(n, 2u);
	unsigned raw = 10 * n;
	Array<double> raw_pf(raw, 2u);

	std::vector< double > xv(30);

	unsigned i = 0;
	for(i = 1; i < 30; i++) xv[i] = 0.;
	i = 0;

	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{
		xv[0] = x;
		raw_pf(i, 0) = ZDT1F1( xv );
		raw_pf(i, 1) = ZDT1F2( xv );
	}
	i = raw - 1;
	xv[0] = xmax;
	raw_pf(i, 0) = ZDT1F1( xv );
	raw_pf(i, 1) = ZDT1F2( xv );

	approxEquidistantFront(raw_pf, pf, n);
}


//************************************************************************
// ZDT 2
//************************************************************************

double ZDT2F1( const std::vector< double >& x )
{
	return x[0];
}

double ZDT2G(const std::vector< double >& x )
{
	double g = 0;
	unsigned n = x.size( );
	for ( unsigned i = 1; i < n; i++ ) 
		g += x[i];
	g = 9 * g / (n - 1.) + 1.;
	return g;
}

double ZDT2F2( const std::vector< double >& x )
{
	return  ZDT2G(x) * (1 - sqr( x[0] /  ZDT2G(x) ));
}

void ZDT2SampleFront(Array< double > &pf, unsigned n = 500) 
{
	double xmin = 0;
	double xmax = 1.;

	pf.resize(n, 2u);
	unsigned raw = 10 * n;
	Array<double> raw_pf(raw, 2u);

	std::vector< double > xv(30);

	unsigned i = 0;
	for(i = 1; i < 30; i++) xv[i] = 0.;
	i = 0;
	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{
		xv[0] = x;
		raw_pf(i, 0) = ZDT2F1( xv );
		raw_pf(i, 1) = ZDT2F2( xv );
	}
	i = raw - 1;
	xv[0] = xmax;
	raw_pf(i, 0) = ZDT2F1( xv );
	raw_pf(i, 1) = ZDT2F2( xv );

	approxEquidistantFront(raw_pf, pf, n);
}

//************************************************************************
// ZDT 3
//************************************************************************

double ZDT3F1( const std::vector< double >& x )
{
	return x[0];
}

double ZDT3G(const std::vector< double >& x )
{
	double g = 0;
	unsigned n = x.size( );
	for ( unsigned i = 1; i < n; i++ ) 
		g += x[i];
	g = 9 * g / (n - 1.) + 1.;
	return g;
}

double ZDT3F2( const std::vector< double >& x )
{
	double gx = ZDT3G(x);
	return  gx * (1 - sqrt( x[0] / gx ) - x[0]/gx * sin (10 * M_PI *  x[0]) );
}

void ZDT3SampleFront(Array< double > &pf, unsigned n = 500) 
{
	double xmin = 0;
	double xmax = 1.;
	std::vector< double > xv(30);

	pf.resize(n, 2u);
	unsigned raw = 40 * n;
	Array<double> raw_pf(raw, 2u);

	unsigned i = 0;
	for(i = 1; i < 30; i++) xv[i] = 0.;
	i = 0;
	double it = (xmax - xmin) / (raw - 1.);
	
	for(double x = xmin; x < (xmax - 0.1 * it); x += it, i++) 
	  {
	    xv[0] = x;
	    raw_pf(i, 0) = ZDT3F1( xv );
	    raw_pf(i, 1) = ZDT3F2( xv );
	  }
	i = raw - 1;
	xv[0] = xmax;
	raw_pf(i, 0) = ZDT3F1( xv );
	raw_pf(i, 1) = ZDT3F2( xv );
	
	sort2DBy1st(raw_pf);
	for(i = 0;  i < (raw_pf.dim(0) - 1); ) {
	  if(raw_pf(i + 1,1) > raw_pf(i, 1)) raw_pf.remove_row(i + 1);
	  else i++;
	}

	approxEquidistantFront(raw_pf, pf, n);
}


//************************************************************************
// ZDT 4
//************************************************************************

// Scaled version
double ZDT4F1( const std::vector< double >& x )
{
  if(x[0] < -5.) {
    std::cerr << "ZDT4F1 called with invalid parameter x[0] = " <<  x[0] << std::endl;
    exit(EXIT_FAILURE);
  }
  return  (x[0] + 5.)/10.;// x[0];
}

double ZDT4G(const std::vector< double >& x )
{
	double g;
	unsigned n = x.size( );
	g = 1 + 10 * (n - 1);
	for ( unsigned i = 1; i < n; i++ ) 
	  g += sqr(x[i]) - 10 * cos(4 * M_PI * x[i]);
	return g;
}

double ZDT4F2( const std::vector< double >& x )
{
  return ZDT4G(x) * ( 1 - sqrt( ((x[0]+5.)/10.)/ ZDT4G(x) ) );
}

double ZDT4PrimeF2( const std::vector< double >& x , const Array<double> &coord)
{
  unsigned i, j, n = x.size();
  double sum;
  std::vector< double > y;
  y.push_back(x[0]);
  for(i = 0; i < n-1; i++) {
    sum = 0.;
    for(j = 1; j < n; j++) sum += coord(i,j - 1) * x[j];
    y.push_back(sum);
  }
  return ZDT4G(y) * ( 1 - sqrt( ((y[0]+5.)/10.)/ ZDT4G(y) ) );
}

// Scaled version 2
double ZDT4FII1( const std::vector< double >& x )
{
  return  (x[0] + 5.)/10.;// x[0];
}

double ZDT4GII(const std::vector< double >& x )
{
	double g;
	unsigned n = x.size( );
	g = 1 + 10 * (n - 1);
	for ( unsigned i = 1; i < n; i++ ) 
	  g += sqr(x[i]) - 10 * cos(4 * M_PI * x[i]);
	return g;
}

double ZDT4FII2( const std::vector< double >& x )
{
  return ZDT4G(x) * ( 1 - sqrt( ((x[0]+5.)/10.)/ ZDT4G(x) ) );
}

// Unscaled version
double ZDT4FG(const std::vector< double >& x )
{
	double g;
	unsigned n = x.size( );
	g = 1 + 10 * (n - 1);
	for ( unsigned i = 1; i < n; i++ ) 
	  g += sqr(x[i]) - 10 * cos(4 * M_PI * x[i]);
	return g;
}

double ZDT4FF2( const std::vector< double >& x )
{
	return  ZDT4FG(x) * (1 - sqrt( x[0] /  ZDT4FG(x) ));
}

double ZDT4FF1( const std::vector< double >& x )
{
	return  x[0];
}


void ZDT4SampleFront(Array< double > &pf, unsigned n = 500) 
{
	double xmin = 0.;
	double xmax = 1.;

	pf.resize(n, 2u);
	unsigned raw = 10 * n;
	Array<double> raw_pf(raw, 2u);

	std::vector< double > xv(10);

	unsigned i = 0;
	for(i = 1; i < 10; i++) xv[i] = 0.;
	i = 0;
	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{
		xv[0] = x;
		raw_pf(i, 0) = ZDT4FF1( xv );
		raw_pf(i, 1) = ZDT4FF2( xv );
	}
	i = raw - 1;
	xv[0] = xmax;
	raw_pf(i, 0) = ZDT4FF1( xv );
	raw_pf(i, 1) = ZDT4FF2( xv );

	approxEquidistantFront(raw_pf, pf, n);
}

//************************************************************************
// ZDT 6
//************************************************************************

double ZDT6F1( const std::vector< double >& x )
{
	return 1-exp(-4*x[0]) * pow( sin(6 * M_PI * x[0]), 6 );
}

double ZDT6G(const std::vector< double >& x )
{
	double g;
	unsigned n = x.size( );
	g = 0;
	for ( unsigned i = 1; i < n; i++ ) 
		g += x[i];
	g /= (n - 1);
	g = pow(g, 0.25);
	g *= 9.;
	g += 1.;
	return g;
}

double ZDT6F2( const std::vector< double >& x )
{
	return  ZDT6G(x) * (1 - (ZDT6F1(x) /  ZDT6G(x)) * (ZDT6F1(x) /  ZDT6G(x)) );
}

void ZDT6SampleFront(Array< double > &pf, unsigned n = 500) 
{
	double xmin = 0.;
	double xmax = 1.;

	pf.resize(n, 2u);
	unsigned raw = 50 * n;
	Array<double> raw_pf(raw, 2u);

	std::vector< double > xv(10);

	unsigned i = 0;
	for(i = 1; i < 10; i++) xv[i] = 0.;
	i = 0;
	double it = (xmax - xmin) / (raw - 1.);

	for(double x = xmin; x < xmax - 0.1 * it; x += it, i++) 
	{
		xv[0] = x;
		raw_pf(i, 0) = ZDT6F1( xv );
		raw_pf(i, 1) = ZDT6F2( xv );
	}
	i = raw - 1;
	xv[0] = xmax;
	raw_pf(i, 0) = ZDT6F1( xv );
	raw_pf(i, 1) = ZDT6F2( xv );

	approxEquidistantFront(raw_pf, pf, n);
}
//
//************************************************************************
// End of ZDT test function implementation of Christian Igel
//************************************************************************
//

////////////////////////////////////////////////////////////////////// SOO
//
//************************************************************************
// SOO test function implementation of Tatsuya Okabe
//************************************************************************
//

//************************************************************************
// Sphere Test Function
//************************************************************************
double sphere( const std::vector< double >& x )
{
	double a;
	unsigned i, n;
	for ( i = 0, a = 0, n = x.size( ); i < n; i++ )
	{
		a += x[ i ]*x[ i ];
	}
	return a;
}

//************************************************************************
// DeJong F2 Test Function
//************************************************************************
double DeJongF2( const std::vector< double >& x )
{
	double f;
	unsigned i, n;
	f = 0.0;
	n = x.size( );
	for ( i = 0; i < n-1; i++ )
	{
		f += 100.0 * pow( x[i]*x[i]-x[i+1], 2.0 ) + pow( 1.0-x[0], 2.0 );
	}
	return f; 
}

//************************************************************************
// DeJong F3 Test Function
//************************************************************************
double DeJongF3( const std::vector< double >& x )
{
	double f;
	unsigned i, n;
	f = 0.0;
	n = x.size( );
	for ( i = 0; i < n; i++ )
	{
		f += (double)( floor( x[ i ] ) );
	}
	return f; 
}

//************************************************************************
// Schaffer F7 Test Function
//************************************************************************
double SchafferF7( const std::vector< double >& x )
{
	double f,t;
	unsigned i, n;
	f = 0.0;
	n = x.size( );
	for ( i = 0; i < n-1; i++ )
	{
		t = x[i]*x[i] + x[i+1]*x[i+1];
		f += pow(t,0.25) * ( pow( sin( 50*pow(t,0.1) ), 2.0 ) + 1.0 );
	}
	return f;
}

//************************************************************************
// Schwefel F1 Test Function
//************************************************************************
double SchwefelF1( const std::vector< double >& x )
{
	double f;
	unsigned i, n;
	f = 0.0;
	n = x.size( );
	for ( i = 0; i < n; i++ )
	{
		f += -x[i] * sin( pow( fabs(x[i]), 0.5 ) );
	}
	return f; 
}

//************************************************************************
// Schwefel F2 Test Function
//************************************************************************
double SchwefelF2( const std::vector< double >& x )
{
	double f;
	unsigned i, j, n;
	f = 0.0;
	n = x.size( );
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < i+1; j++ )
		{
			f += x[j] * x[j];
		}
	}
	return f; 
}

//************************************************************************
// Rastrigin Test Function
//************************************************************************
double rastrigin( const std::vector< double >& x )
{
	unsigned i, n;
	double sum;
	const double C = znPi2C;
	double A = 10.; // suitable for dim = 10
	double B = 10.; // suitable for dim = 10
	for ( i = 0, n = x.size( ), sum = n*A; i < n; i++ )
	{
		sum += ( x[ i ] * x[ i ] ) - B * cos( C * x[ i ] );
	}
	return sum;
}

//************************************************************************
// Rosenbrock Test Function
//************************************************************************
double rosenbrock( const std::vector< double >& x )
{
	const double A = 100.;
	const double B = 1.;
	unsigned i, n;
	double   a;
	for( i = 0, a = 0, n = x.size( ); i < n; i++ )
	{
		a += A*( ( ( ( x[ i ]) - x[ i+1 ]) * ( (x[ i ]) - x[ i+1 ]) ) + ( (x[ i ] - B) * (x[ i ] - B) ) ) * ( ( ( ( x[ i ]) - x[ i+1 ]) * ( (x[ i ]) - x[ i+1 ]) ) + ( (x[ i ] - B) * (x[ i ] - B) ) );
	}  
	return a;
}

//************************************************************************
// Ackley Test Function
//************************************************************************
double ackley( const std::vector< double >& x )
{
	const double A        = 20.;
	const double B        = 0.2;
	const double C        = znPi2C;
	unsigned i, n;
	double   a, b;
	for ( a = b = 0., i = 0, n = x.size( ); i < n; ++i )
	{
		a += x[ i ] * x[ i ];
		b += cos( C * x[ i ] );
	}
	return -A * exp( -B * sqrt( a / n ) ) - exp( b / n ) + A + znEC;
}

//************************************************************************
// Func1 Test Function ( no name )
//************************************************************************
double func1( const std::vector< double >& x )
{
	const double A = 5.;
	const double B = 31.4159265359;
	unsigned i, n;
	double   a;
	for ( i = 0, a  = 0, n = x.size( ); i < n; i++ )
	{
		if ( (x[ i ]>=-5) && (x[ i ]<=5) )
		{
			a += (A - fabs(x [ i ])) * fabs( cos( B * x [ i ]) );
		}
	}
	return a;
}




//
//************************************************************************
// IHR
//************************************************************************
//

//////////
// IHR1 //
//////////
inline double h(double x, double n) {
  return 1/(1+exp(-x/sqrt(n)));
}

inline double hf(double x, double y0, double ymax) {
  if(fabs(y0) <=  ymax) return x;
  return fabs(y0) + 1;
}

inline double hg(double x) {
  return sqr(x) / (fabs(x) + 0.1);
}

double IHR1F1(const std::vector<double> &x, const Array<double> &coord)
{
  return fabs(scalarprod(coord.row(0), x));
}

double IHR1F2(const std::vector<double> &x, const Array<double> &coord)
{
  unsigned i; 
  double g = 0;
  unsigned n = x.size();
  Array<double> y(n);
  
  for(i = 0; i < n; i++) y(i) = scalarprod(coord.row(i), x);
  
  double ymax = fabs(coord.col(0)(0));
  for(i = 1; i < n; i++) if(fabs(coord.row(0)(i)) > ymax) ymax = fabs(coord.row(0)(i));
  ymax = 1/ymax;
  
  for (i = 1; i < n; i++ ) g += hg(y(i));
  g = 9 * g / (n - 1.) + 1.;

  if(fabs(y(0)) <= ymax) return g * (1. - sqrt(h(y(0), n) / g));
  return g * (1 + fabs(y(0)));
}

//////////
// IHR2 //
//////////
double IHR2F1(const std::vector<double> &x, Array<double> &coord)
{
	return fabs(scalarprod(coord.col(0), x));
}

double IHR2F2(const std::vector<double> &x, Array<double> &coord)
{
  unsigned i;
  double g = 0;
  unsigned n = x.size();
  Array<double> y;
  y.resize(n);
  
  for(i = 0; i < n; i++) y(i) = scalarprod(coord.col(i), x);  
  double ymax = fabs(coord.col(0)(0));
  for(i = 1; i < n; i++) if(fabs(coord.row(0)(i)) > ymax) ymax = fabs(coord.row(0)(i));
  ymax = 1/ymax;
  
  
  for (i = 1; i < n; i++ ) g += hg(y(i));
  g = 9 * g / (n - 1.) + 1.;
  
  if(fabs(y(0)) <= ymax) return g * (1. - sqr(y(0) / g));
  return g * (1 + fabs(y(0)));
}

//////////
// IHR3 //
//////////
double IHR3F1(const std::vector<double> &x, Array<double> &coord)
{
	return fabs(scalarprod(coord.col(0), x));
}

double IHR3F2(const std::vector<double> &x, Array<double> &coord)
{
  unsigned i;
  double g = 0;
  unsigned n = x.size();
  Array<double> y;
  y.resize(n);
  
  for(i = 0; i < n; i++) y(i) = scalarprod(coord.col(i), x);  
  double ymax = fabs(coord.col(0)(0));
  for(i = 1; i < n; i++) if(fabs(coord.row(0)(i)) > ymax) ymax = fabs(coord.row(0)(i));
  ymax = 1/ymax;
  
  
  for (i = 1; i < n; i++ ) g += hg(y(i));
  g = 9 * g / (n - 1.) + 1.;
  
  if(fabs(y(0)) <= ymax) return g * (1. - sqrt(h(y(0), n) / g) - h(y(0), n) / g * sin(10 * M_PI * y(0)));
  return g * (1 + fabs(y(0)));
}


//////////
// IHR4 //
//////////
double IHR4F1(const std::vector<double> &x, Array<double> &coord)
{
  return fabs(scalarprod(coord.col(0), x));
}

double IHR4F2(const std::vector<double> &x, Array<double> &coord)
{
  unsigned i; 
  double g = 0;
  unsigned n = x.size();
  Array<double> y;
  y.resize(n);
  
  for(i = 0; i < n; i++) y(i) = scalarprod(coord.row(i), x);
  
  double ymax = fabs(coord.col(0)(0));
  for(i = 1; i < n; i++) if(fabs(coord.row(0)(i)) > ymax) ymax = fabs(coord.row(0)(i));
  ymax = 1/ymax;
  
  for (i = 1; i < n; i++ ) g += sqr(y(i)) - 10 * cos(4*M_PI*y(i));
  g += 10 * (n - 1.) + 1.;

  if(fabs(y(0)) <= ymax) return g * (1. - sqrt(h(y(0), n) / g));
  return g * (1 + fabs(y(0)));
}


//////////
// IHR6 //
//////////
double IHR6F1(const std::vector<double> &x, Array<double> &coord)
{
  double y0 = scalarprod(coord.col(0), x);
  return 1-exp(-4*fabs(y0)) * pow( sin(6 * M_PI * y0), 6 );
}

double IHR6F2(const std::vector<double> &x, Array<double> &coord)
{
  unsigned i; 
  double g = 0;
  unsigned n = x.size();
  Array<double> y;
  y.resize(n);
  
  for(i = 0; i < n; i++) y(i) = scalarprod(coord.row(i), x);
  double ymax = fabs(coord.col(0)(0));
  for(i = 1; i < n; i++) if(fabs(coord.row(0)(i)) > ymax) ymax = fabs(coord.row(0)(i));
  ymax = 1/ymax;
    
  for (i = 1; i < n; i++ ) g += hg(y(i));
  g = 1 + 9 * pow(g / (n-1.), .25);

  if(fabs(y(0)) <= ymax) return g * (1. - sqr(IHR6F1(x, coord) / g));
  return g * (1 + fabs(y(0)));
}


#endif /* !__TESTFUNCTIOMOO_H */

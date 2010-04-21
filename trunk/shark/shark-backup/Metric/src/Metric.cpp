/*! ======================================================================
 *
 *  File    :  Metric.cpp
 *  Created :  2000-09-05
 *
 *  Copyright (c) 2000 Axel W. Dietrich
 *
 *  \par Copyright (c) 1999-2003:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>

 *  \par Project:
 *      Metric
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: Metric.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: Metric.cpp,v $
 *      Revision 2.1  2004/06/17 23:52:59  saviapbe
 *      The doxygen's documentation generation for Metric was added.
 *
 *      Revision 2.0  2003/11/28 16:23:11  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *<BR>
 
 * ----------------------------------------------------------------------
 *
 *  This file is part of the Metric. This library is free software;
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

#include "Metric/Metric.h"

/*!
* \par Compute the Levenshtein distance of two vectors.
* 
*
* 
*/
int 
Metric::Levenshtein( const vector<bool> &a, const vector<bool> &b )
{
  register int   i, j;
  int            dL;
  int            rows= a.size();
  int            columns= b.size();
  int          **table= new int*[ rows ];

  for ( i= 0; i<rows; i++ ) {
	table[ i ] = new int[ columns ];
	table[ i ][ 0 ]= i;
  }
  for ( j= 0; j<columns; j++ )
	table[ 0 ][ j ]= j;

  for( j= 1; j<columns; j++ ) {
	for( i= 1; i<rows; i++ ) {
	  table[ i ][ j ]= min( 1 + table[i][j-1], 1 + table[i-1][j],
						    !(a[i-1]==b[j-1]) + table[i-1][j-1]
						  );
	}
  }
  dL= table[ rows-1 ][ columns-1 ];

  for ( i= 0; i<rows; i++ )
	delete [ ] table[ i ];
  delete [] table;

  return( dL );
}


/*!
* Compute the Levenshtein distance of two (char*) bitstrings. 
* Make sure that the bitstrings have the form "00110101010...".
*
* \example metric_test.cpp
*/
int
Metric::Levenshtein( const char *a, const char *b )
{
  register unsigned i;
  vector< bool > va( strlen( a ) );
  vector< bool > vb( strlen( b ) );
  
  // copy da strings into the vectors
  for ( i= 0; i<va.size(); ++i )
	va[ i ]= a[ i ]-'0';
  for ( i= 0; i<vb.size(); ++i )
	vb[ i ]= b[ i ]-'0';
  
  return Levenshtein( va, vb );
}


//
// Compute the Levenshtein distance of two strings. Strings 
// must look like "01001001010...".
//
static int Levenshtein( const string &a, const string &b )
{
  const char *ca= const_cast<char*>( a.c_str() );
  const char *cb= const_cast<char*>( b.c_str() );
  
  return Levenshtein( ca, cb );
}


// 
// Compute the Levenshtein distance of two integer arrays.
// Arrays must contain a -1 as the last element, eg. { 0,1,1,0,-1 }.
//
int
Metric::Levenshtein( const int *a, const int *b )
{
  register unsigned i;
  vector< bool > va( arrayLength( a ) + 1 );
  vector< bool > vb( arrayLength( b ) + 1 );
  
  // copy da integer Arrays into the vectors
  for ( i= 0; i<va.size(); ++i )
	va[ i ]= a[ i ];
  for ( i= 0; i<vb.size(); ++i )
	vb[ i ]= b[ i ];
  
  return Levenshtein( va, vb );
}


// ======================================================================


/*!
* \par Compute the Haming distance of two (boolean, bitstring) vectors.
*
* \example metric_test.cpp
*/
int 
Metric::Hamming( const vector<bool>& a, const vector<bool>& b )
{
  int  x= a.size();
  int  y= b.size();
  int dH= 0;
  
  if ( x!=y || x==0 || y==0) {
	  cerr << "Vectors are of different size or one/both is/are mt. Can't compute Hamming distance!\n";
	  return -1;
  }
  while( x-- )
	if ( a[--y] != b[y] ) 
	  dH++;
  return( dH );
}


//
// Compute the Hamming distance dH of two strings. Make
// sure that the strings look like "00110101010...".
//
int
Metric::Hamming( const string &a, const string &b )
{
  char *ca= const_cast<char*>( a.c_str() );
  char *cb= const_cast<char*>( b.c_str() );
  
  return Hamming( ca, cb );
}


/*!
* \par Compute the Hamming distance dH of two (char *) bitstrings. Make
* sure that the bitstrings look like "00110101010...".
* 
* \example metric_test.cpp
*/
int 
Metric::Hamming( const char *a, const char *b )
{
  register unsigned i;
  vector< bool > va( strlen( a ) );
  vector< bool > vb( strlen( b ) );
  
  // copy da strings into the vectors
  for ( i= 0; i<va.size(); ++i )
	va[ i ]= a[ i ]-'0';
  for ( i= 0; i<vb.size(); ++i )
	vb[ i ]= b[ i ]-'0';
  
  return Hamming( va, vb );
}


//
// Compute the Hamming distance dH of two integer arrays.
// Arrays must contain a -1 as the last element, eg. { 0,1,1,0,-1 }.
//
int
Metric::Hamming( const int *a, const int *b )
{
  register unsigned i;
  vector< bool > va( arrayLength( a ) + 1 );
  vector< bool > vb( arrayLength( b ) + 1 );
  
  // copy da strings into the vectors
  for ( i= 0; i<va.size(); ++i )
	va[ i ]= a[ i ];
  for ( i= 0; i<vb.size(); ++i )
	vb[ i ]= b[ i ];
  
  return Hamming( va, vb );
}



// ======================================================================

//
// Compute the minimum of three integers.
//
int
Metric::min( const int a, const int b, const int c ) 
{
  if ( a<=b  &&  a<=c ) return a;
  if ( b<=a  &&  b<=c ) return b;
  return c;
}


//
// Compute the number of elements in an integer array. End-of-array
// marker must be a -1.
//
int
Metric::arrayLength( const int *array ) 
{
  register int i;		
  for ( i= 0; *array!=-1; ++i, ++array ) 
	/* mt */ ;
  return i;
}





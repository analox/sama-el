/*! Individual.cpp
* ======================================================================
*
*  File    :  Individual.cpp
*  Created :  01.01.1995
*
*  Copyright (c) 1995-2000 Martin Kreutz
*
*  Institut fuer Neuroinformatik
*  Ruhr-Universitaet Bochum
*  44780 Bochum, Germany<BR>
*        Phone: +49-234-32-25558<BR>
*        Fax:   +49-234-32-14209<BR>
*        eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
*        www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
*  
*    \par Project:
*        EALib
*  
*    \par Language and Compiler:
*        C++, egcs (Linux)
*  
*    \par File and Revision:
*        $RCSfile: Individual.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: Individual.cpp,v $
*        Revision 2.6  2004/06/17 23:09:00  saviapbe
*        An error in the INT's URL was corrected.
*
*        Revision 2.5  2004/06/17 22:53:19  saviapbe
*        The standard file header was for doxygen adapted.
*
*        Revision 2.4  2004/05/20 12:41:26  shark-admin
*        pvm_pkind() and pvm_upkind() modified
*  
*        Revision 2.3  2004/03/05 16:39:24  shark-admin
*        change of some member names: initObjBuffer(...) => initDisplVal(...), getNumObj() => getNoOfDisplVal(), setObjective(...) => setDisplVal(...), getObjective(...) => getDisplVal(...), objective => displVal, numObj => noDisplVal.
*  
*        Revision 2.2  2004/01/30 15:24:31  shark-admin
*        PVM interface modified
*  
*        Revision 2.1  2004/01/12 11:06:01  shark-admin
*        buffer for objectives (indepently from fitness values) modified
*        pvm pack routines modified
*  
*        Revision 2.0  2003/11/28 16:23:09  shark-admin
*        Revision tag reset to revision tag 2.x
*  
*        Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
*        INI Administration
*  <BR>
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

#include "EALib/Individual.h"

using namespace std;

//===========================================================================
//
// constructors
//
Individual::Individual( unsigned n )
  : vector< Chromosome * >( n )
{
    for( unsigned i = size( ); i--; )
        *( begin( ) + i ) = new ChromosomeT< char >;

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( unsigned n, const Chromosome& chrom )
  : vector< Chromosome * >( n )
{
    for( unsigned i = size( ); i--; )
        *( begin( ) + i ) = chrom.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;    
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0 )
  : vector< Chromosome * >( 1 )
{
    *( begin( ) ) = chrom0.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0,
		        const Chromosome& chrom1 )
  : vector< Chromosome * >( 2 )
{
    *( begin( )     ) = chrom0.clone( );
    *( begin( ) + 1 ) = chrom1.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0,
		        const Chromosome& chrom1,
		        const Chromosome& chrom2 )
  : vector< Chromosome * >( 3 )
{
    *( begin( )     ) = chrom0.clone( );
    *( begin( ) + 1 ) = chrom1.clone( );
    *( begin( ) + 2 ) = chrom2.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0,
		        const Chromosome& chrom1,
		        const Chromosome& chrom2,
		        const Chromosome& chrom3 )
  : vector< Chromosome * >( 4 )
{
    *( begin( )     ) = chrom0.clone( );
    *( begin( ) + 1 ) = chrom1.clone( );
    *( begin( ) + 2 ) = chrom2.clone( );
    *( begin( ) + 3 ) = chrom3.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0,
		        const Chromosome& chrom1,
		        const Chromosome& chrom2,
		        const Chromosome& chrom3,
		        const Chromosome& chrom4 )
  : vector< Chromosome * >( 5 )
{
    *( begin( )     ) = chrom0.clone( );
    *( begin( ) + 1 ) = chrom1.clone( );
    *( begin( ) + 2 ) = chrom2.clone( );
    *( begin( ) + 3 ) = chrom3.clone( );
    *( begin( ) + 4 ) = chrom4.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0,
		        const Chromosome& chrom1,
		        const Chromosome& chrom2,
		        const Chromosome& chrom3,
		        const Chromosome& chrom4,
		        const Chromosome& chrom5 )
  : vector< Chromosome * >( 6 )
{
    *( begin( )     ) = chrom0.clone( );
    *( begin( ) + 1 ) = chrom1.clone( );
    *( begin( ) + 2 ) = chrom2.clone( );
    *( begin( ) + 3 ) = chrom3.clone( );
    *( begin( ) + 4 ) = chrom4.clone( );
    *( begin( ) + 5 ) = chrom5.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0,
		        const Chromosome& chrom1,
		        const Chromosome& chrom2,
		        const Chromosome& chrom3,
		        const Chromosome& chrom4,
		        const Chromosome& chrom5,
		        const Chromosome& chrom6 )
  : vector< Chromosome * >( 7 )
{
    *( begin( )     ) = chrom0.clone( );
    *( begin( ) + 1 ) = chrom1.clone( );
    *( begin( ) + 2 ) = chrom2.clone( );
    *( begin( ) + 3 ) = chrom3.clone( );
    *( begin( ) + 4 ) = chrom4.clone( );
    *( begin( ) + 5 ) = chrom5.clone( );
    *( begin( ) + 6 ) = chrom6.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Chromosome& chrom0,
		        const Chromosome& chrom1,
		        const Chromosome& chrom2,
		        const Chromosome& chrom3,
		        const Chromosome& chrom4,
		        const Chromosome& chrom5,
		        const Chromosome& chrom6,
		        const Chromosome& chrom7 )
  : vector< Chromosome * >( 8 )
{
    *( begin( )     ) = chrom0.clone( );
    *( begin( ) + 1 ) = chrom1.clone( );
    *( begin( ) + 2 ) = chrom2.clone( );
    *( begin( ) + 3 ) = chrom3.clone( );
    *( begin( ) + 4 ) = chrom4.clone( );
    *( begin( ) + 5 ) = chrom5.clone( );
    *( begin( ) + 6 ) = chrom6.clone( );
    *( begin( ) + 7 ) = chrom7.clone( );

    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const vector< Chromosome * >& chrom )
  : vector< Chromosome * >( chrom )
{
    fitness       = 0;
    scaledFitness = 0;
    evalFlg       = true;
    feasible      = 0;
    selProb       = 0.;
    numCopies     = 0;
    elitist       = false;
    age           = 0;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = 0;
    noDisplVal    = 0;
    //
}

Individual::Individual( const Individual& indiv )
  : vector< Chromosome * >( indiv.size( ) )
{
    for( unsigned i = size( ); i--; )
        *( begin( ) + i ) = indiv[ i ].clone( );

    fitness       = indiv.fitness;
    scaledFitness = indiv.scaledFitness;
    evalFlg       = indiv.evalFlg;
    feasible      = indiv.feasible;
    selProb       = indiv.selProb;
    numCopies     = indiv.numCopies;
    elitist       = indiv.elitist;
    age           = indiv.age;
    //added by Stefan Wiegand 24.02.2002
    learnTime     = indiv.learnTime;
    noDisplVal    = indiv.getNoOfDisplVal();
    //added by Stefan Wiegand 12.01.2004
    Individual::initDisplVal(noDisplVal);
    for (unsigned ii = 0; ii<noDisplVal; ii++ )
      setDisplVal(ii, indiv.getDisplVal(ii));
    //
}

//===========================================================================
//
// destructor
//
Individual::~Individual ()
{
    for( unsigned i = size( ); i--; )
        delete *( begin( ) + i );
}

//===========================================================================

unsigned Individual::totalSize( ) const
{
    unsigned s = 0;
    for( unsigned i = size( ); i--; s += ( *( begin( ) + i ) )->size( ) );
    return s;
}

//===========================================================================
// changed by Marc Toussaint and Stefan Wiegand (INI) 20.11.2002
Individual& Individual::operator = ( const Individual& indiv )
{
    unsigned i;

    for( i = size( ); i--; )
        delete *( begin( ) + i );
    vector< Chromosome * >::operator = ( indiv );
    for( i = size( ); i--; )
        *( begin( ) + i ) = indiv[ i ].clone( );
  
    // added by Marc Toussaint and Stefan Wiegand (INI) 20.11.2002
    for( i = size( ); i--; )
      (*( begin( ) + i ))->registerIndividual(*this,i);
    //

    fitness       = indiv.fitness;
    scaledFitness = indiv.scaledFitness;
    evalFlg       = indiv.evalFlg;
    feasible      = indiv.feasible;
    selProb       = indiv.selProb;
    numCopies     = indiv.numCopies;
    elitist       = indiv.elitist;
    age           = indiv.age;
    learnTime     = indiv.learnTime;
    noDisplVal    = indiv.getNoOfDisplVal();
 
    // added by Stefan Wiegand (INI) 12.01.2004
    Individual::initDisplVal(noDisplVal);
    for (i = 0; i<noDisplVal; i++ )
      setDisplVal(i, indiv.getDisplVal(i));
    //

    return *this;
}

//===========================================================================

void Individual::replace( unsigned i, const Chromosome& chrom )
{
    RANGE_CHECK( i < size( ) )
    delete *( begin( ) + i );
    *( begin( ) + i ) = chrom.clone( );
}

void Individual::insert( unsigned i, const Chromosome& chrom )
{
    RANGE_CHECK( i <= size( ) )
    vector< Chromosome * >::insert( begin( ) + i, chrom.clone( ) );
}

void Individual::append( const Chromosome& chrom )
{
    vector< Chromosome * >::push_back( chrom.clone( ) );
}

void Individual::remove( unsigned i )
{
    RANGE_CHECK( i < size( ) )
    delete *( begin( ) + i );
    vector< Chromosome * >::erase( begin( ) + i );
}

void Individual::remove( unsigned i, unsigned k )
{
    if( i <= k ) {
        RANGE_CHECK( k < size( ) )
	for( unsigned j = i; j < k; j++ )
	    delete *( begin( ) + j );
	vector< Chromosome * >::erase( begin( ) + i, begin( ) + k );
    }
}

//===========================================================================

bool Individual::operator == ( const Individual& ind ) const
{
    if( size( ) == ind.size( ) ) {
        for( unsigned i = 0; i < size( ); ++i )
	    if( ! ( ( *this )[ i ] == ind[ i ] ) ) return false;

	return fitness       == ind.fitness       &&
	  scaledFitness == ind.scaledFitness &&
	  evalFlg       == ind.evalFlg       &&
	  feasible      == ind.feasible      &&
	  selProb       == ind.selProb       &&
	  numCopies     == ind.numCopies     &&
	  elitist       == ind.elitist       &&
	  age           == ind.age           &&
	  learnTime     == ind.learnTime    ;
    }
    return false;
}

bool Individual::operator < ( const Individual& ind ) const
{
    if( size( ) == ind.size( ) ) {
        bool less = false;

        for( unsigned i = 0; i < size( ); ++i )
	    if( ind[ i ] < ( *this )[ i ] )
	        return false;
	    else if( ! less && ( *this )[ i ] < ind[ i ] )
	        less = true;

	return less/*                             &&
	       fitness       <= ind.fitness       &&
	       scaledFitness <= ind.scaledFitness &&
	       evalFlg       <= ind.evalFlg       &&
	       feasible      <= ind.feasible      &&
	       selProb       <= ind.selProb       &&
	       numCopies     <= ind.numCopies     &&
	       elitist       <= ind.elitist       &&
	       age           <= ind.age*/;
    }

    return size( ) < ind.size( );
}

//===========================================================================
// added by Marc Toussaint and Stefan Wiegand (INI) 20.11.2002 

int Individual::pvm_pkind(){
  //cout << "\t Individual_pk" << endl;

  unsigned i;

  unsigned *s = new unsigned;
  *s = this->size();
  pvm_pkuint(s,1,1);
  delete s; 

  for (i = 0; i < this->size(); i++)
    ((*this)[i]).pvm_pkchrom();

  unsigned *u = new unsigned[7];
  u[0] = feasible;
  u[1] = elitist;
  u[2] = evalFlg;
  u[3] = numCopies;
  u[4] = age;
  u[5] = learnTime;
  u[6] = getNoOfDisplVal();
  pvm_pkuint(u,7,1);
  delete[] u; 

  double *f = new double[displVal.size()+3];
  unsigned jj=0;
  for (;jj<displVal.size();jj++)
  f[jj]   = getDisplVal(jj); 
  f[jj]   = fitness;
  f[++jj] = scaledFitness;
  f[++jj] = selProb;
  pvm_pkdouble(f,displVal.size()+3,1);
  delete[] f; 
  
  return 1;
}

int Individual::pvm_upkind(){
  //cout << "\t Individual_upk" << endl;

  unsigned i;

  unsigned *s = new unsigned;
  pvm_upkuint(s,1,1);
  if(this->size()!=*s){
    std::cerr << "EALib/Individual.cpp: the individual which has called pvm_upkind() is of unexpected size!" << std::endl;
    std::cerr << "Please initialize this individual prototypically with chromosomes of appropriate type." << std::endl;
    exit(-1);
  }
  delete s; 

  for (i = 0; i < this->size(); i++)
    ((*this)[i]).pvm_upkchrom();
  
  unsigned *u = new unsigned[7];
  pvm_upkuint(u,7,1);
  feasible  = (u[0] != 0);
  elitist   = (u[1] != 0);
  evalFlg   = (u[2] != 0);
  numCopies = u[3];
  age       = u[4];
  learnTime = u[5];
  noDisplVal= u[6];
  delete[] u;

  double *f = new double[noDisplVal+3];
  pvm_upkdouble(f,noDisplVal+3,1);
  unsigned jj=0;
  for (;jj< noDisplVal ;jj++)
    setDisplVal(jj,f[jj]); 
  fitness       = f[jj];
  scaledFitness = f[++jj];
  selProb       = f[++jj];
  delete[] f;
  
  return 1;
}

//===========================================================================

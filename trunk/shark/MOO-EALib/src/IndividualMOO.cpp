/* ======================================================================
 *
 *  File    :  IndividualMOO.cpp
 *  Created :  2002-01-23
 *
 *  Description : Class IndividualMOO in EALib for MOO
 *  Author  : Tatsuya Okabe <tatsuya.okabe@honda-ri.de>
 *  Copyright (c) 2002-2004 GPL2 Honda Research Institute Europe GmbH
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
 *      $RCSfile: IndividualMOO.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: IndividualMOO.cpp,v $
 *      Revision 2.12  2005/01/31 12:20:21  shark-admin
 *      declarations removed from for(...) header
 *
 *      Revision 2.11  2005/01/31 09:42:48  shark-admin
 *      definitions of UnpenalizedFitness functions added
 *
 *      Revision 2.10  2004/11/15 15:37:06  shark-admin
 *      compatibility with vc++ ensured
 *
 *      Revision 2.9  2004/06/18 03:45:50  saviapbe
 *      Dummy linking to examples.
 *
 *      Revision 2.8  2004/05/20 12:05:08  shark-admin
 *      pvm_pkind( ) and pvm_upkind( ) modified
 *
 *      Revision 2.7  2004/04/27 15:06:51  shark-admin
 *      *** empty log message ***
 *
 *      Revision 2.5  2004/04/12 08:17:05  shark-admin
 *       IndividualMOO::size( ) and IndividualMOO::totalSize( ) modified
 *
 *      Revision 2.4  2004/03/05 16:43:30  shark-admin
 *      change of some member names in Individual.h: initDisplVal(...), getNoOfDisplVal(), setDisplVal(...), getDisplVal(...).
 *
 *      Revision 2.3  2004/03/03 13:23:49  saviapbe
 *      Copyright information was changed.
 *
 *      Revision 2.2  2004/02/12 11:27:09  saviapbe
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

#include "MOO-EALib/IndividualMOO.h"
#include "EALib/Individual.h"
#include <typeinfo>

using namespace std;

// ----------------------------------------------------------------------------
// TO-IM-001
/*!
*
*
*
* \example VEGA.cpp
*
*/
IndividualMOO::IndividualMOO( )
  : Individual( )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}  
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-002
IndividualMOO::IndividualMOO( unsigned n )
     : Individual( n )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-003
IndividualMOO::IndividualMOO( unsigned n, const Chromosome& chrom1 ) 
     : Individual( n, chrom1 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-004
IndividualMOO::IndividualMOO( const Chromosome& chrom1 ) 
     : Individual( chrom1 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-005
IndividualMOO::IndividualMOO( const Chromosome& chrom1,
                              const Chromosome& chrom2 ) 
     : Individual( chrom1, chrom2 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-006
IndividualMOO::IndividualMOO( const Chromosome& chrom1,
                              const Chromosome& chrom2,
                              const Chromosome& chrom3 ) 
     : Individual( chrom1, chrom2, chrom3 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-007
IndividualMOO::IndividualMOO( const Chromosome& chrom1,
                              const Chromosome& chrom2,
                              const Chromosome& chrom3,
                              const Chromosome& chrom4 ) 
     : Individual( chrom1, chrom2, chrom3, chrom4 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-008
IndividualMOO::IndividualMOO( const Chromosome& chrom1,
                              const Chromosome& chrom2,
                              const Chromosome& chrom3,
                              const Chromosome& chrom4,
                              const Chromosome& chrom5 ) 
     : Individual( chrom1, chrom2, chrom3, chrom4, chrom5 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-009
IndividualMOO::IndividualMOO( const Chromosome& chrom1,
                              const Chromosome& chrom2,
                              const Chromosome& chrom3,
                              const Chromosome& chrom4,
                              const Chromosome& chrom5,
                              const Chromosome& chrom6 ) 
     : Individual( chrom1, chrom2, chrom3, chrom4, chrom5, chrom6 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-010
IndividualMOO::IndividualMOO( const Chromosome& chrom1,
                              const Chromosome& chrom2,
                              const Chromosome& chrom3,
                              const Chromosome& chrom4,
                              const Chromosome& chrom5,
                              const Chromosome& chrom6,
                              const Chromosome& chrom7 ) 
     : Individual( chrom1, chrom2, chrom3, chrom4, chrom5, chrom6, chrom7 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-011
IndividualMOO::IndividualMOO( const Chromosome& chrom1,
                              const Chromosome& chrom2,
                              const Chromosome& chrom3,
                              const Chromosome& chrom4,
                              const Chromosome& chrom5,
                              const Chromosome& chrom6,
                              const Chromosome& chrom7,
                              const Chromosome& chrom8 ) 
     : Individual( chrom1, chrom2, chrom3, chrom4, chrom5, chrom6, chrom7, chrom8 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-012
IndividualMOO::IndividualMOO( const vector< Chromosome* >& chrom1 )
     : Individual( chrom1 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-013
IndividualMOO::IndividualMOO( const Individual& indiv1 ) 
     : Individual( indiv1 )
{
  MOOFitness.resize( 1 );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
  setMOORank( 0 );
  setMOOShare( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-014
IndividualMOO::IndividualMOO( const IndividualMOO& indmoo ) 
  : Individual( dynamic_cast< const Individual& >( indmoo ) )
{
  unsigned i;
  // copy data
  MOOFitness.resize( indmoo.MOOFitness.size( ) );
  for ( i = indmoo.MOOFitness.size( ); i--; ){
    setMOOFitness( i, indmoo.MOOFitness[i] );
  }
  UnpenalizedMOOFitness.resize( indmoo.UnpenalizedMOOFitness.size( ) );
  for ( i = indmoo.UnpenalizedMOOFitness.size( ); i--; ){
    setUnpenalizedMOOFitness( i, indmoo.UnpenalizedMOOFitness[i] );
  }
  setMOORank( indmoo.MOORank );
  setMOOShare( indmoo.MOOShare );
  fitness       = indmoo.fitness;
  scaledFitness = indmoo.scaledFitness;
  evalFlg       = indmoo.evalFlg;
  feasible      = indmoo.feasible;
  selProb       = indmoo.selProb;
  numCopies     = indmoo.numCopies;
  elitist       = indmoo.elitist;
  age           = indmoo.age;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-015
IndividualMOO::~IndividualMOO( ){ }

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-020
unsigned IndividualMOO::size( ) const
{
  //return vector< Chromosome * >::size( );
  return  Individual::size( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-021
unsigned IndividualMOO::totalSize( ) const
{
  //   unsigned s = 0;

  //   for ( unsigned i = size( ); i--; ){
  //     s += (*(vector< Chromosome *>::begin( ) + i ))->size( );
  //   }
  //   return s;
  return Individual::totalSize();
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-030
void IndividualMOO::setFitness( double fit )
{
  Individual::setFitness( fit );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-031
double IndividualMOO::fitnessValue( ) const
{
  return Individual::fitnessValue( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-032
double IndividualMOO::getFitness( ) const
{
  return Individual::fitnessValue( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-033
void IndividualMOO::setScaledFitness( double scalef )
{
  Individual::setScaledFitness( scalef );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-034
double IndividualMOO::getScaledFitness( ) const
{
  return Individual::getScaledFitness( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-035
void IndividualMOO::setAge( unsigned age )
{
  Individual::setAge( age );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-036
void IndividualMOO::incAge( )
{
  Individual::setAge( Individual::getAge( ) + 1 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-037
unsigned IndividualMOO::getAge( ) const
{
  return Individual::getAge( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-038
void IndividualMOO::setSelectionProbability( double prob )
{
  Individual::setSelectionProbability( prob );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-039
void IndividualMOO::setSelProb( double prob )
{
  Individual::setSelectionProbability( prob );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-040
double IndividualMOO::selectionProbability( ) const
{
  return Individual::selectionProbability( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-041
double IndividualMOO::getSelProb( ) const
{
  return Individual::selectionProbability( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-042
void IndividualMOO::setNumCopies( unsigned num )
{ 
  Individual::setNumCopies( num );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-043
unsigned IndividualMOO::numberOfCopies( ) const
{
  return Individual::numberOfCopies( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-044
unsigned IndividualMOO::getNumCopies( ) const
{
  return Individual::numberOfCopies( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-045
void IndividualMOO::setEvaluationFlag( )
{
  Individual::setEvaluationFlag( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-046
void IndividualMOO::clearEvaluationFlag( )
{
  Individual::clearEvaluationFlag( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-047
void IndividualMOO::setEvalFlg( bool flg )
{
  if ( flg ){
    Individual::setEvaluationFlag( );
  }
  else{
    Individual::clearEvaluationFlag( );
  }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-048
bool IndividualMOO::needEvaluation( ) const
{
  return Individual::needEvaluation( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-049
bool IndividualMOO::getEvalFlg( ) const
{
  return Individual::needEvaluation( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-050
void IndividualMOO::setFeasible( bool fea )
{
  Individual::setFeasible( fea );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-051
bool IndividualMOO::isFeasible( ) const
{
  return Individual::isFeasible( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-052
bool IndividualMOO::getFeasible( ) const
{
  return Individual::isFeasible( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-053
void IndividualMOO::setElitist( bool eli )
{
  Individual::setElitist( eli );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-054
bool IndividualMOO::isElitist( ) const
{
  return Individual::isElitist( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-055
bool IndividualMOO::getElitist( ) const
{
  return Individual::isElitist( );
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// SW-IM-050
void IndividualMOO::setLearnTime( unsigned lt )
{
  Individual::setLearnTime( lt );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SW-IM-051
unsigned IndividualMOO::getLearnTime( ) const
{
  return Individual::getLearnTime( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SW-IM-052
unsigned IndividualMOO::getNoOfDisplVal() const 
{ 
  return Individual::getNoOfDisplVal();
};
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SW-IM-053
void IndividualMOO::initDisplVal(unsigned o)
{
  Individual::initDisplVal(o);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SW-IM-054
void IndividualMOO::setDisplVal( unsigned i, double e)
{ 
  Individual::setDisplVal( i , e );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SW-IM-055
double IndividualMOO::getDisplVal( unsigned i ) const
{
 return Individual::getDisplVal( i );
} 
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-060
void IndividualMOO::setNoOfObj( unsigned n )
{
  MOOFitness.resize( n );
  UnpenalizedMOOFitness.resize( n );
  initializeMOOFitness( 0.0 );
  initializeUnpenalizedMOOFitness( 0.0 );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-061
unsigned IndividualMOO::getNoOfObj( )
{
  return MOOFitness.size( );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-062
void IndividualMOO::setMOORank( unsigned n )
{
  MOORank = n;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-063
unsigned IndividualMOO::getMOORank( )
{
  return MOORank;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-064
void IndividualMOO::setMOOShare( double n )
{
  MOOShare = n;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-065
double IndividualMOO::getMOOShare( )
{
  return MOOShare;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-066
void IndividualMOO::setMOOFitness( unsigned nof, double fit )
{
  RANGE_CHECK( nof < MOOFitness.size( ) ) 
  MOOFitness[nof] = fit;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-067
double IndividualMOO::getMOOFitness( unsigned nof )
{
  RANGE_CHECK( nof < MOOFitness.size( ) )
  return MOOFitness[nof];
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-068
void IndividualMOO::setMOOFitnessValues( double f0 )
{
  RANGE_CHECK( 0 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-069
void IndividualMOO::setMOOFitnessValues( double f0, double f1 )
{
  RANGE_CHECK( 1 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
  MOOFitness[1] = f1;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-070
void IndividualMOO::setMOOFitnessValues( double f0, double f1,
                                         double f2 )
{
  RANGE_CHECK( 2 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
  MOOFitness[1] = f1;
  MOOFitness[2] = f2;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-071
void IndividualMOO::setMOOFitnessValues( double f0, double f1,
                                         double f2, double f3 )
{
  RANGE_CHECK( 3 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
  MOOFitness[1] = f1;
  MOOFitness[2] = f2;
  MOOFitness[3] = f3;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-072
void IndividualMOO::setMOOFitnessValues( double f0, double f1,
                                         double f2, double f3,
                                         double f4 )
{
  RANGE_CHECK( 4 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
  MOOFitness[1] = f1;
  MOOFitness[2] = f2;
  MOOFitness[3] = f3;
  MOOFitness[4] = f4;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-073
void IndividualMOO::setMOOFitnessValues( double f0, double f1,
                                         double f2, double f3,
                                         double f4, double f5 )
{
  RANGE_CHECK( 5 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
  MOOFitness[1] = f1;
  MOOFitness[2] = f2;
  MOOFitness[3] = f3;
  MOOFitness[4] = f4;
  MOOFitness[5] = f5;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-074
void IndividualMOO::setMOOFitnessValues( double f0, double f1,
                                         double f2, double f3,
                                         double f4, double f5,
                                         double f6 )
{
  RANGE_CHECK( 6 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
  MOOFitness[1] = f1;
  MOOFitness[2] = f2;
  MOOFitness[3] = f3;
  MOOFitness[4] = f4;
  MOOFitness[5] = f5;
  MOOFitness[6] = f6;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-075
void IndividualMOO::setMOOFitnessValues( double f0, double f1,
                                         double f2, double f3,
                                         double f4, double f5,
                                         double f6, double f7 )
{
  RANGE_CHECK( 7 < MOOFitness.size( ) )
  MOOFitness[0] = f0;
  MOOFitness[1] = f1;
  MOOFitness[2] = f2;
  MOOFitness[3] = f3;
  MOOFitness[4] = f4;
  MOOFitness[5] = f5;
  MOOFitness[6] = f6;
  MOOFitness[7] = f7;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-076
void IndividualMOO::setMOOFitnessValues( vector< double >& fit )
{
  RANGE_CHECK( fit.size( ) < MOOFitness.size( ) )
  for ( unsigned i = fit.size( ); i--; ){
    MOOFitness[i] = fit[i];
  }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-077
vector< double >& IndividualMOO::getMOOFitnessValues( )
{
  return MOOFitness;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-078
void IndividualMOO::initializeMOOFitness( double x )
{
  for ( unsigned i = MOOFitness.size( ); i--; ){
    MOOFitness[i] = x;
  }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-079
void IndividualMOO::setUnpenalizedMOOFitness( unsigned nof, double fit )
{
  RANGE_CHECK( nof < UnpenalizedMOOFitness.size( ) ) 
  UnpenalizedMOOFitness[nof] = fit;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-080
double IndividualMOO::getUnpenalizedMOOFitness( unsigned nof )
{
  RANGE_CHECK( nof < UnpenalizedMOOFitness.size( ) )
  return UnpenalizedMOOFitness[nof];
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-081
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0 )
{
  RANGE_CHECK( 0 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-082
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0, double f1 )
{
  RANGE_CHECK( 1 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
  UnpenalizedMOOFitness[1] = f1;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-083
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0, double f1,
						  double f2 )
{
  RANGE_CHECK( 2 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
  UnpenalizedMOOFitness[1] = f1;
  UnpenalizedMOOFitness[2] = f2;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-084
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0, double f1,
						  double f2, double f3 )
{
  RANGE_CHECK( 3 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
  UnpenalizedMOOFitness[1] = f1;
  UnpenalizedMOOFitness[2] = f2;
  UnpenalizedMOOFitness[3] = f3;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-085
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0, double f1,
						  double f2, double f3,
						  double f4 )
{
  RANGE_CHECK( 4 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
  UnpenalizedMOOFitness[1] = f1;
  UnpenalizedMOOFitness[2] = f2;
  UnpenalizedMOOFitness[3] = f3;
  UnpenalizedMOOFitness[4] = f4;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-086
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0, double f1,
						  double f2, double f3,
						  double f4, double f5 )
{
  RANGE_CHECK( 5 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
  UnpenalizedMOOFitness[1] = f1;
  UnpenalizedMOOFitness[2] = f2;
  UnpenalizedMOOFitness[3] = f3;
  UnpenalizedMOOFitness[4] = f4;
  UnpenalizedMOOFitness[5] = f5;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-087
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0, double f1,
						  double f2, double f3,
						  double f4, double f5,
						  double f6 )
{
  RANGE_CHECK( 6 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
  UnpenalizedMOOFitness[1] = f1;
  UnpenalizedMOOFitness[2] = f2;
  UnpenalizedMOOFitness[3] = f3;
  UnpenalizedMOOFitness[4] = f4;
  UnpenalizedMOOFitness[5] = f5;
  UnpenalizedMOOFitness[6] = f6;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-088
void IndividualMOO::setUnpenalizedMOOFitnessValues( double f0, double f1,
						  double f2, double f3,
						  double f4, double f5,
						  double f6, double f7 )
{
  RANGE_CHECK( 7 < UnpenalizedMOOFitness.size( ) )
  UnpenalizedMOOFitness[0] = f0;
  UnpenalizedMOOFitness[1] = f1;
  UnpenalizedMOOFitness[2] = f2;
  UnpenalizedMOOFitness[3] = f3;
  UnpenalizedMOOFitness[4] = f4;
  UnpenalizedMOOFitness[5] = f5;
  UnpenalizedMOOFitness[6] = f6;
  UnpenalizedMOOFitness[7] = f7;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-089
void IndividualMOO::setUnpenalizedMOOFitnessValues( vector< double >& fit )
{
  RANGE_CHECK( fit.size( ) < UnpenalizedMOOFitness.size( ) )
  for ( unsigned i = fit.size( ); i--; ){
    UnpenalizedMOOFitness[i] = fit[i];
  }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-090
vector< double >& IndividualMOO::getUnpenalizedMOOFitnessValues( )
{
  return UnpenalizedMOOFitness;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SR-IM-091
void IndividualMOO::initializeUnpenalizedMOOFitness( double x )
{
  for ( unsigned i = UnpenalizedMOOFitness.size( ); i--; ){
    UnpenalizedMOOFitness[i] = x;
  }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-100
Chromosome& IndividualMOO::operator [ ] ( unsigned i )
{
  return Individual::operator [ ] ( i );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-101
const Chromosome& IndividualMOO::operator [ ] ( unsigned i ) const
{
  return Individual::operator [ ] ( i );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-102
IndividualMOO& IndividualMOO::operator = ( const IndividualMOO& indmoo )
{
  //
  // replaced by Stefan Wiegand (INI) 14.1.2004
  //
  unsigned i;
  
  for( i = size( ); i--; )
    delete *( begin( ) + i );

  vector< Chromosome * >::operator = ( indmoo );

  for( i = size( ); i--; )
    *( begin( ) + i ) = indmoo[ i ].clone( );
  
  for( i = size( ); i--; )
    (*( begin( ) + i ))->registerIndividual(*this,i);
  //
  //   if ( size( ) != indmoo.size( ) ){
  //     cout << "\n ***** The size is unmatched in TO-IM-102 *****\n" << endl;
  //     cout << size( ) << " <--- " << indmoo.size( ) << endl;
  //     cout << "\n Please set the same number of Chromosome *****\n" << endl;
  //     return *this;
  //   }
  //
  if ( MOOFitness.size( ) != indmoo.MOOFitness.size( ) ){
    setNoOfObj( indmoo.MOOFitness.size( ) );
  }
  
  //for ( unsigned i = size( ); i--; ){
  //  // *(*( vector< Chromosome * >::begin( ) + i )) = indmoo[i];
  //  *(*( begin( ) + i )) = indmoo[i];
  //}
 
  (*this).fitness       = indmoo.fitness;
  (*this).scaledFitness = indmoo.scaledFitness;
  (*this).evalFlg       = indmoo.evalFlg;
  (*this).feasible      = indmoo.feasible;
  (*this).selProb       = indmoo.selProb;
  (*this).numCopies     = indmoo.numCopies;
  (*this).elitist       = indmoo.elitist;
  (*this).age           = indmoo.age;
  (*this).learnTime     = indmoo.getLearnTime();
  (*this).noDisplVal    = indmoo.getNoOfDisplVal();
  (*this).MOORank       = indmoo.MOORank;
  (*this).MOOShare      = indmoo.MOOShare;

  for ( i=getNoOfObj( ); i--; ){
    MOOFitness[i] = indmoo.MOOFitness[i];
  }

  for ( i=getNoOfObj( ); i--; ){
    UnpenalizedMOOFitness[i] = indmoo.UnpenalizedMOOFitness[i];
  }

  //
  // added by Stefan Wiegand (INI) 12.01.2004
  Individual::initDisplVal( Individual::noDisplVal);
  for (i = 0; i< Individual::noDisplVal; i++ ){
    setDisplVal(i, indmoo.getDisplVal(i));
  }
  return *this;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-103
bool IndividualMOO::operator == ( const IndividualMOO& indmoo ) const
{
  unsigned i;
  if ( size( ) == indmoo.size( ) ){
    // ***** check Chromosome
    for ( i = size( ); i--; ){
      if ( !( (*this)[i] == indmoo[i] ) ){
	return false;
      }
    }
    // ***** check MOO fitness
    if ( MOOFitness.size( ) != indmoo.MOOFitness.size( ) || UnpenalizedMOOFitness.size( )!= indmoo.UnpenalizedMOOFitness.size( )){
      return false;
    }

    for ( i = MOOFitness.size( ); i--; ){
      if ( !( (*this).MOOFitness[i] == indmoo.MOOFitness[i] ) ){
	return false;
      }
    }

    for ( i = UnpenalizedMOOFitness.size( ); i--; ){
      if ( !( (*this).UnpenalizedMOOFitness[i] == indmoo.UnpenalizedMOOFitness[i] ) ){
	return false;
      }
    }

    // ***** check internal variables
    if ( fitness       == indmoo.fitness       &&
         scaledFitness == indmoo.scaledFitness &&
         evalFlg       == indmoo.evalFlg       &&
         feasible      == indmoo.feasible      &&
         selProb       == indmoo.selProb       &&
         numCopies     == indmoo.numCopies     &&
         elitist       == indmoo.elitist       &&
         age           == indmoo.age           &&
         MOORank       == indmoo.MOORank       &&
         MOOShare      == indmoo.MOOShare         ){
      return true;
    }
  }
  return false;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-104
IndividualMOO& IndividualMOO::operator = ( const Individual& ind )
{
  //
  // replaced by Stefan Wiegand (INI) 14.1.2004
  //
  unsigned i;
  
  for( i = size( ); i--; )
    delete *( begin( ) + i );

  vector< Chromosome * >::operator = ( ind );

  for( i = size( ); i--; )
    *( begin( ) + i ) = ind[ i ].clone( );
  
  for( i = size( ); i--; )
    (*( begin( ) + i ))->registerIndividual(*this,i);
  //
  //   if ( size( ) != ind.size( ) ){
  //     cout << "\n ***** The size is unmatched in TO-IM-104 *****\n" << endl;
  //     cout << size( ) << " <--- " << ind.size( ) << endl;
  //     cout << "\n Please set the same number of Chromosome *****\n" << endl;
  //     return *this;
  //   }
  //
  setNoOfObj( 1 );
  initializeMOOFitness( 0 );
  initializeUnpenalizedMOOFitness( 0 );

  //
  //   for ( unsigned i = size( ); i--; ){
  //   //     *(*( vector< Chromosome * >::begin( ) + i )) = ind[i];
  //     *(*( begin( ) + i )) = ind[i];
  //   }

  (*this).fitness       = ind.fitnessValue( );
  (*this).scaledFitness = ind.getScaledFitness( );
  (*this).evalFlg       = ind.needEvaluation( );
  (*this).feasible      = ind.isFeasible( );
  (*this).selProb       = ind.selectionProbability( );
  (*this).numCopies     = ind.numberOfCopies( );
  (*this).elitist       = ind.isElitist( );
  (*this).age           = ind.getAge( );    
  (*this).learnTime     = ind.getLearnTime( );
  (*this).noDisplVal    = ind.getNoOfDisplVal();
  (*this).MOORank       = 0;
  (*this).MOOShare      = 0;

  // added by Stefan Wiegand (INI) 12.01.2004
  Individual::initDisplVal(noDisplVal);
  for (i = 0; i<noDisplVal; i++ )
    setDisplVal(i, ind.getDisplVal(i));
  
  return *this;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-110
void IndividualMOO::replace( unsigned i, const Chromosome& chrom )
{
  Individual::replace( i, chrom );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-111
void IndividualMOO::insert( unsigned i, const Chromosome& chrom )
{
  Individual::insert( i, chrom );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-112
void IndividualMOO::append( const Chromosome& chrom )
{
  //Individual::append( chrom ); 
  vector< Chromosome * >::push_back( chrom.clone( ) );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-113
void IndividualMOO::remove( unsigned i )
{
  Individual::remove( i );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-114
void IndividualMOO::remove( unsigned i, unsigned j )
{
  Individual::remove( i, j );
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-200
double IndividualMOO::aggregation( const vector< double >& weight )
{
  if ( (*this).MOOFitness.size( ) > weight.size( ) ){
    cout << "\n ***** The size is not match in TO-IM-200 *****\n" << endl;
    return 0.0;
  }
  double sum = 0.0;
  for ( unsigned i = (*this).MOOFitness.size( ); i--; ){
    sum += (*this).getMOOFitness( i ) * weight[ i ];
  }
  (*this).setFitness( sum );
  return sum;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-201
double IndividualMOO::simplesum( )
{
  double sum = 0.0;
  for ( unsigned i = (*this).MOOFitness.size( ); i--; ){
    sum += (*this).getMOOFitness( i );
  }
  (*this).setFitness( sum );
  return sum;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TO-IM-500
void IndividualMOO::printIM( )
{
  unsigned i;
  cout << "\n\n********** IndividualMOO **********\n";
  cout << "No. of Chromosomes         : " << (*this).size( ) << "\n";
  cout << "Variable ( fitness )       : " << (*this).fitness << "\n";
  cout << "Variable ( scaledFitness ) : " << (*this).scaledFitness << "\n";
  cout << "Variable ( evalFlg )       : " << (*this).evalFlg << "\n";
  cout << "Variable ( feasible )      : " << (*this).feasible << "\n";
  cout << "Variable ( selProb )       : " << (*this).selProb << "\n";
  cout << "Variable ( numCopies )     : " << (*this).numCopies << "\n";
  cout << "Variable ( elitist )       : " << (*this).elitist << "\n";
  cout << "Variable ( age )           : " << (*this).age << "\n";
  cout << "Variable ( MOOShare )      : " << (*this).MOOShare << "\n";
  cout << "Variable ( MOORank )       : " << (*this).MOORank << "\n";
  cout << "No. of fitness functions   : " << (*this).MOOFitness.size( ) 
                                          << "\n";
  for ( i = 0; i < (*this).MOOFitness.size( ); i++ ){
    cout << "Value of fitness function  : " << (*this).MOOFitness[i]
         << " ( fun. = " << i << " )\n";
  }

  for ( i = 0; i < (*this).UnpenalizedMOOFitness.size( ); i++ ){
    cout << "Value of penalized fitness function  : " << (*this).UnpenalizedMOOFitness[i]
         << " ( fun. = " << i << " )\n";
  }

  if ( (*this).size( ) != 0 ){
    for ( unsigned i = 0; i < size( ); i++ ){
      cout << "No. Of Alleles             : " << (*this)[i].size( ) 
           << " ( Chr. = " << i << " )\n";
    }
  }
  cout << endl;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SW-IM-600
int IndividualMOO::pvm_pkind()
{
  //cout << "\t IndividualMOO_pk" << endl;

  unsigned i;

  unsigned *s = new unsigned;
  *s = this->size();
  pvm_pkuint(s,1,1);
  delete s; 

  for (i = 0; i < this->size(); i++)
    ((*this)[i]).pvm_pkchrom();
  
  uint *u = new uint[9];
  u[0] = feasible;
  u[1] = elitist;
  u[2] = evalFlg;
  u[3] = numCopies;
  u[4] = age;
  u[5] = learnTime;
  u[6] = getNoOfDisplVal();
  u[7] = MOORank;
  u[8] = getNoOfObj();
  pvm_pkuint(u,9,1);
  delete[] u; 

  double *g = new double[2*getNoOfObj()];
  unsigned jj=0;
  for (;jj<getNoOfObj();jj++)
    g[jj]=MOOFitness[jj];
  for (;jj<2*getNoOfObj();jj++)
    g[jj]=UnpenalizedMOOFitness[jj];
  pvm_pkdouble(g,2*getNoOfObj(),1);
  delete[] g; 
  
  double *f = new double[getNoOfDisplVal()+4];
  jj=0;
  for (;jj<getNoOfDisplVal();jj++)
  f[jj]   = getDisplVal(jj); 
  f[jj]   = fitness;
  f[++jj] = scaledFitness;
  f[++jj] = selProb;
  f[++jj] = MOOShare;
  pvm_pkdouble(f,getNoOfDisplVal()+4,1);
  delete[] f; 
  
  return 1;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// SW-IM-601
int IndividualMOO::pvm_upkind()
{
  //cout << "\t Individual_upk" << endl;

  unsigned i;

  unsigned *s = new unsigned;
  pvm_upkuint(s,1,1);
  if(this->size()!=*s){
    std::cerr << "MOO-EALib/IndividualMOO.cpp: the individual which has called pvm_upkind() is of unexpected size!" << std::endl;
    std::cerr << "Please initialize this individual prototypically with chromosomes of appropriate type." << std::endl;
    exit(-1);
  }
  delete s; 

  for (i = 0; i < this->size(); i++)
    ((*this)[i]).pvm_upkchrom();
  
  uint *o = new uint;
  uint *u = new uint[9];
  pvm_upkuint(u,9,1);
  feasible  = u[0]?1:0;
  elitist   = u[1]?1:0;
  evalFlg   = u[2]?1:0;
  numCopies = u[3];
  age       = u[4];
  learnTime = u[5];
  noDisplVal= u[6];
  MOORank   = u[7];
  *o         = u[8];
  delete[] u;

  setNoOfObj(*o);

  double *g = new double[2 * getNoOfObj()];
  pvm_upkdouble(g,2 * getNoOfObj(),1); 
  unsigned jj=0;
  for (;jj<getNoOfObj();jj++)
    MOOFitness[jj]=g[jj];
  for (;jj<2 * getNoOfObj();jj++)
    UnpenalizedMOOFitness[jj]=g[jj];
  delete[] g; 
  delete[] o;

  double *f = new double[getNoOfDisplVal()+4];
  pvm_upkdouble(f,getNoOfDisplVal()+4,1);
  jj=0;
  for (;jj< getNoOfDisplVal() ;jj++)
    setDisplVal(jj,f[jj]); 
  fitness       = f[jj];
  scaledFitness = f[++jj];
  selProb       = f[++jj];
  MOOShare      = f[++jj];
  delete[] f;

  return 1;
}
// ----------------------------------------------------------------------------






/* ======================================================================
 *
 *  File    :  IndividualMOO.h
 *  Created :  2001-12-16
 *
 *  Description : Class IndividualMOO in EALib for MOO
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
 *      $RCSfile: IndividualMOO.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: IndividualMOO.h,v $
 *      Revision 2.8  2005/06/10 13:03:31  glasmtbl
 *      *** empty log message ***
 *
 *      Revision 2.7  2005/01/14 15:03:04  shark-admin
 *      additional objective vector 'UnpenalizedMOOFitness' added
 *
 *      Revision 2.6  2004/11/15 15:36:19  shark-admin
 *      compatibility with vc++ ensured
 *
 *      Revision 2.5  2004/05/20 12:06:49  shark-admin
 *      documentation modified
 *
 *      Revision 2.4  2004/03/05 16:44:09  shark-admin
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



#ifndef __INDIVIDUALMOO_H
#define __INDIVIDUALMOO_H


#include "EALib/Individual.h"
#include "EALib/ChromosomeFactory.h"


class IndividualMOO : public Individual
{

 public:

  // ---------------------------------------------------  
  // constructor
  // ---------------------------------------------------
  // TO-IM-001
  IndividualMOO( );
  // TO-IM-002
  IndividualMOO( unsigned );
  // TO-IM-003
  IndividualMOO( unsigned, const Chromosome& );
  // TO-IM-004
  IndividualMOO( const Chromosome& );
  // TO-IM-005
  IndividualMOO( const Chromosome&,
		 const Chromosome& );
  // TO-IM-006
  IndividualMOO( const Chromosome&,
		 const Chromosome&,
		 const Chromosome& );
  // TO-IM-007
  IndividualMOO( const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome& );
  // TO-IM-008
  IndividualMOO( const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome& );
  // TO-IM-009
  IndividualMOO( const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome& );
  // TO-IM-010
  IndividualMOO( const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome& );
  // TO-IM-011
  IndividualMOO( const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome&,
		 const Chromosome& );
  // TO-IM-012
  IndividualMOO( const std::vector< Chromosome* >& );
  // TO-IM-013
  IndividualMOO( const Individual& );
  // TO-IM-014
  IndividualMOO( const IndividualMOO& );
  
  // ---------------------------------------------------  
  // destructor
  // ---------------------------------------------------
  // TO-IM-015
  ~IndividualMOO( );

  // ---------------------------------------------------  
  // structure
  // ---------------------------------------------------
  // TO-IM-020
  unsigned size( ) const;
  // TO-IM-021
  unsigned totalSize( ) const;

  // ---------------------------------------------------  
  // internal variables
  // ---------------------------------------------------

  // fitness
  // TO-IM-030
  void setFitness( double );
  // TO-IM-031
  double fitnessValue( ) const;
  // TO-IM-032
  double getFitness( ) const;

  // scaledFitness
  // TO-IM-033
  void setScaledFitness( double );
  // TO-IM-034
  double getScaledFitness( ) const;

  // age
  // TO-IM-035
  void setAge( unsigned );
  // TO-IM-036
  void incAge( );
  // TO-IM-037
  unsigned getAge( ) const;

  // selProb
  // TO-IM-038
  void setSelectionProbability( double );
  // TO-IM-039
  void setSelProb( double );
  // TO-IM-040
  double selectionProbability( ) const;
  // TO-IM-041
  double getSelProb( ) const;

  // numCopies
  // TO-IM-042
  void setNumCopies( unsigned );
  // TO-IM-043
  unsigned numberOfCopies( ) const;
  // TO-IM-044
  unsigned getNumCopies( ) const;
  
  // evalFlag
  // TO-IM-045
  void setEvaluationFlag( );
  // TO-IM-046
  void clearEvaluationFlag( );
  // TO-IM-047
  void setEvalFlg( bool );
  // TO-IM-048
  bool needEvaluation( ) const;
  // TO-IM-049
  bool getEvalFlg( ) const;

  // feasible
  // TO-IM-050
  void setFeasible( bool );
  // TO-IM-051
  bool isFeasible( ) const;
  // TO-IM-052
  bool getFeasible( ) const;

  // elitist
  // TO-IM-053
  void setElitist( bool );
  // TO-IM-054
  bool isElitist( ) const;
  // TO-IM-055
  bool getElitist( ) const;

  // parameter for embedded learning
  // SW-IM-050
  void setLearnTime( unsigned lt = 0 );
  // SW-IM-051
  unsigned getLearnTime( ) const; 

  //Interface for a buffer that stores performance values
  //indepently from fitness
  // SW-IM-052
  unsigned getNoOfDisplVal( ) const;
  // SW-IM-053
  void initDisplVal( unsigned o=0 );
  // SW-IM-054
  void setDisplVal( unsigned i, double e = 0 );
  // SW-IM-055
  double getDisplVal( unsigned i ) const;

  // ---------------------------------------------------  
  // internal variables ( IndividualMOO )
  // ---------------------------------------------------

  // Number of objectives
  // TO-IM-060
  void setNoOfObj( unsigned );
  // TO-IM-061
  unsigned getNoOfObj( );

  // MOORank
  // TO-IM-062
  void setMOORank( unsigned );
  // TO-IM-063
  unsigned getMOORank( );

  // MOOShare
  // TO-IM-064
  void setMOOShare( double );
  // TO-IM-065
  double getMOOShare( );

  // MOOFitness
  // TO-IM-066
  void setMOOFitness( unsigned, double );
  // TO-IM-067
  double getMOOFitness( unsigned );
  // TO-IM-068
  void setMOOFitnessValues( double );
  // TO-IM-069
  void setMOOFitnessValues( double, double );
  // TO-IM-070
  void setMOOFitnessValues( double, double, double );
  // TO-IM-071
  void setMOOFitnessValues( double, double, double,
			    double );
  // TO-IM-072
  void setMOOFitnessValues( double, double, double,
			    double, double );
  // TO-IM-073
  void setMOOFitnessValues( double, double, double,
			    double, double, double );
  // TO-IM-074
  void setMOOFitnessValues( double, double, double,
			    double, double, double,
			    double );
  // TO-IM-075
  void setMOOFitnessValues( double, double, double,
			    double, double, double,
			    double, double );
  // TO-IM-076
  void setMOOFitnessValues( std::vector< double >& );
  // TO-IM-077
  std::vector< double >& getMOOFitnessValues( );
  // TO-IM-078
  void initializeMOOFitness( double = 0.0 );


  // Unpenalized MOOFitness
  // SR-IM-079
  void setUnpenalizedMOOFitness( unsigned, double );
  // SR-IM-080
  double getUnpenalizedMOOFitness( unsigned );
  // SR-IM-081
  void setUnpenalizedMOOFitnessValues( double );
  // SR-IM-082
  void setUnpenalizedMOOFitnessValues( double, double );
  // SR-IM-083
  void setUnpenalizedMOOFitnessValues( double, double, double );
  // SR-IM-084
  void setUnpenalizedMOOFitnessValues( double, double, double,
			    double );
  // SR-IM-085
  void setUnpenalizedMOOFitnessValues( double, double, double,
			    double, double );
  // SR-IM-086
  void setUnpenalizedMOOFitnessValues( double, double, double,
			    double, double, double );
  // SR-IM-087
  void setUnpenalizedMOOFitnessValues( double, double, double,
			    double, double, double,
			    double );
  // SR-IM-088
  void setUnpenalizedMOOFitnessValues( double, double, double,
			    double, double, double,
			    double, double );
  // SR-IM-089
  void setUnpenalizedMOOFitnessValues( std::vector< double >& );
  // SR-IM-090
  std::vector< double >& getUnpenalizedMOOFitnessValues( );
  // SR-IM-091
  void initializeUnpenalizedMOOFitness( double = 0.0 );
  
  // ---------------------------------------------------  
  // operator
  // ---------------------------------------------------
  // TO-IM-100
  Chromosome& operator [ ] ( unsigned );
  // TO-IM-101
  const Chromosome& operator [ ] ( unsigned ) const;
  // TO-IM-102
  IndividualMOO& operator = ( const IndividualMOO& );
  // TO-IM-103
  bool operator == ( const IndividualMOO& ) const;
  // TO-IM-104
  IndividualMOO& operator = ( const Individual& );

  // ---------------------------------------------------  
  // Change Chromosomes
  // ---------------------------------------------------
  // TO-IM-110
  void replace( unsigned, const Chromosome& );
  // TO-IM-111
  void insert( unsigned, const Chromosome& );
  // TO-IM-112
  void append( const Chromosome& );
  // TO-IM-113
  void remove( unsigned );
  // TO-IM-114
  void remove( unsigned, unsigned );

  // SW-IM-110
  // possibly conflicting with 'to-im-112': see vc-compiler message 
  template<class ChromosomeTemplate> void append( const ChromosomeTemplate& chrom );

  // ---------------------------------------------------  
  // calculation for fitness
  // ---------------------------------------------------
  // TO-IM-200
  double aggregation( const std::vector< double >& );
  // TO-IM-201
  double simplesum( );
  
  // ---------------------------------------------------  
  // Print data
  // ---------------------------------------------------
  // TO-IM-500
  void printIM( );

  // ---------------------------------------------------  
  // PVM Interface
  // ---------------------------------------------------

  // SW-IM-600
  /*! Part of PVM-send routine for MOO-individuals */
  int pvm_pkind();

  // SW-IM-601
  /*! Part of PVM-receive routine for MOO-individuals */
  int pvm_upkind();


#ifndef __NO_GENERIC_IOSTREAM
	friend std::ostream& operator << ( std::ostream& os, const IndividualMOO& ind )
	{
//		printf("IndividualMOO::operator << \n"); fflush(stdout);

		os << "IndividualMOO(" << ind.size( ) << ")\n"
			<< ind.fitness       << '\n'
			<< ind.scaledFitness << '\n'
			<< ind.evalFlg       << '\n'
			<< ind.feasible      << '\n'
			<< ind.selProb       << '\n'
			<< ind.numCopies     << '\n'
			<< ind.elitist       << '\n'
			<< ind.age           << '\n'      
			<< ind.learnTime     << '\n'    
			<< ind.getNoOfDisplVal()	<< '\n'
			<< ind.MOORank				<< '\n'
			<< ind.MOOShare;

		int r, rc;

		rc = ind.MOOFitness.size();
		os << '\n' << rc;
		for (r=0; r<rc; r++)
			os << '\n' << ind.MOOFitness[r];

		rc = ind.UnpenalizedMOOFitness.size();
		os << '\n' << rc;
		for (r=0; r<rc; r++)
			os << '\n' << ind.UnpenalizedMOOFitness[r];

		// added by Stefan Wiegand 07.01.2004
		for(unsigned ii=0; ii<ind.getNoOfDisplVal();ii++)
			os << '\n' << ind.getDisplVal(ii);
		os << std::endl;

		for( unsigned i = 0; i < ind.size( ); ++i )
			os << '\n' << ind[ i ];
		os << std::endl;

		return os;
    }

	friend std::istream& operator >> ( std::istream& is, IndividualMOO& ind )
	{
//		printf("IndividualMOO::operator >> \n"); fflush(stdout);

		unsigned i, indSize( 0 );
		std::string s, t;
		double e;

		is >> s;
		is.get( );   // skip end of line

		if( is.good( ) &&
			s.substr( 0, 14 ) == "IndividualMOO(" &&
			s.find( ')' ) != std::string::npos )
		{
			// Extract the size indication from the string:
			t = s.substr( s.find( '(' )+1, s.find( ')' ) - s.find( '(' ) - 1 );
			indSize = atoi( t.c_str( ) );

			// Adapt size of Individual:
			ind.resize( indSize );

			is >> ind.fitness
				>> ind.scaledFitness
				>> ind.evalFlg
				>> ind.feasible
				>> ind.selProb
				>> ind.numCopies
				>> ind.elitist
				>> ind.age
				>> ind.learnTime
				>> ind.noDisplVal
				>> ind.MOORank
				>> ind.MOOShare;

			int r, rc;

			is >> rc;
			ind.MOOFitness.resize(rc);
			for (r=0; r<rc; r++)
				is >> ind.MOOFitness[r];

			is >> rc;
			ind.UnpenalizedMOOFitness.resize(rc);
			for (r=0; r<rc; r++)
				is >> ind.UnpenalizedMOOFitness[r];

			// added by Stefan Wiegand 07.01.2004
			ind.initDisplVal(ind.getNoOfDisplVal());
			for(i=0; i<ind.getNoOfDisplVal();i++)
			{
				is >> e;
				ind.setDisplVal(i,e);
			}

			for( i = 0; i < ind.size( ); ++i )
			{
				int pos = is.tellg();
				is >> s;
				is.seekg(pos, std::ios_base::beg);
				t = s.substr(s.find('<' )+1, s.find('>') - s.find('<' ) - 1);
				const char* type = t.c_str();

				Chromosome* pC = CreateChromosome(type);
				ind.replace(i, *pC);
				is >> ind[ i ];	
			}
		}
		return is;
	}
#endif // !__NO_GENERIC_IOSTREAM


 protected:
  std::vector< double > MOOFitness;
  std::vector< double > UnpenalizedMOOFitness;
  unsigned         MOORank;
  double           MOOShare;
};




#endif /* !__INDIVIDUALMOO_H */

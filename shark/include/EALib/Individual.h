//===========================================================================
/*!
 *  \file Individual.h
 *
 *
 *  \author  Martin Kreutz
 *  \date    01.01.1995
 *
 *  \par Copyright (c) 1995-2003:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>
 *
 *  \par Project:
 *      EALib
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: Individual.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: Individual.h,v $
 *      Revision 2.5  2004/05/20 12:32:37  shark-admin
 *      documentation enhanced
 *
 *      Revision 2.4  2004/04/12 07:42:19  shark-admin
 *      qualifier of class declaration 'private Individual' changed to 'protected Individual'
 *
 *      Revision 2.3  2004/03/05 16:38:46  shark-admin
 *      change of some member names: initObjBuffer(...) => initDisplVal(...), getNumObj() => getNoOfDisplVal(), setObjective(...) => setDisplVal(...), getObjective(...) => getDisplVal(...), objective => displVal, numObj => noDisplVal.
 *
 *      Revision 2.2  2004/01/30 15:26:02  shark-admin
 *      friend class IndividualMOO added
 *
 *      Revision 2.1  2004/01/12 11:07:03  shark-admin
 *      buffer for objectives (indepently from fitness) modified;
 *      pvm pack routines modified;
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of EALib. This library is free software;
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
 *
 */
//===========================================================================




#ifndef __INDIVIDUAL_H
#define __INDIVIDUAL_H

#include "ChromosomeT.h"
#include "PVMinterface.h"

//===========================================================================

class Individual : protected std::vector< Chromosome * >
{
  public:
    Individual( ) { }
    explicit Individual( unsigned );
    Individual( unsigned, const Chromosome& );
    Individual( const Chromosome& );
    Individual( const Chromosome&,
	        const Chromosome& );
    Individual( const Chromosome&,
	        const Chromosome&,
	        const Chromosome& );
    Individual( const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome& );
    Individual( const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome& );
    Individual( const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome& );
    Individual( const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome& );
    Individual( const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome&,
	        const Chromosome& );
    Individual( const std::vector< Chromosome* >& );
    Individual( const Individual& );
    virtual ~Individual( );

    unsigned size( ) const
	{
		return static_cast< const std::vector< Chromosome * > * >( this )->size( );
	}

    unsigned totalSize           ( ) const;
    double   fitnessValue        ( ) const { return fitness;   }
    double   selectionProbability( ) const { return selProb;   }
    bool     isFeasible          ( ) const { return feasible;  }
    unsigned numberOfCopies      ( ) const { return numCopies; }
    bool     isElitist           ( ) const { return elitist;   }

    void     setFitness             ( double fit ) { fitness=scaledFitness=fit; }
    void     setFeasible            ( bool f ) { feasible = f; }
    void     setSelectionProbability( double ps ) { selProb = ps; }
 
    Chromosome& operator [ ] ( unsigned i )
    {
        RANGE_CHECK( i < size( ) )
		return *( *( begin( ) + i ) );
    }

    const Chromosome& operator [ ] ( unsigned i ) const
    {
        RANGE_CHECK( i < size( ) )
		return *( *( begin( ) + i ) );
    }

    //=======================================================================

    //
    // Changed by Marc Toussaint & Stefan Wiegand at 20.11.2002
    // refer to Individual& Individual::operator = ( const Individual& );
    // in Individual.cpp

    /*! The Chromosome method 'registerIndividual(Ind&,int)' was added, 
      refer also to 'Chromosome.h'.*/

    Individual& operator = ( const Individual& );

    //=======================================================================

    void replace( unsigned i, const Chromosome& chrom );
    void insert ( unsigned i, const Chromosome& chrom );
    void append ( const Chromosome& chrom );
    void remove ( unsigned i );
    void remove ( unsigned from, unsigned to );

    void setAge ( unsigned a = 0 ) { age = a;    }
    void incAge ( )                { ++age;      }
    unsigned getAge ( ) const      { return age; }

    //=======================================================================

    //  Added by Marc Toussaint & Stefan Wiegand at 20.11.2002
    /*! Appends Chromosomes to the Individual, refer also to Chromosome.h. */
  
    template<class ChromosomeTemplate> void append( const ChromosomeTemplate& chrom ){
      ChromosomeTemplate* newChrom=new ChromosomeTemplate(chrom);
      std::vector< Chromosome * >::push_back( newChrom );
      newChrom->registerIndividual(*this,size()-1);
    }


    // Added by Stefan Wiegand at 20.11.2002
    /*! Interface for a buffer that stores how many learning iterations 
        an individual conducts per generation */

    void setLearnTime     ( unsigned lt = 0 ){ learnTime = lt;   }
    unsigned getLearnTime ( ) const          { return learnTime; } 
   

    // Added by Stefan Wiegand at 24.02.2003
    /*! Interface for a flag that indicates the neccessity of an individual to become evaluated */

    void setEvaluationFlag  ( )       { evalFlg = true;  }
    void clearEvaluationFlag( )       { evalFlg = false; }
    bool needEvaluation     ( ) const { return evalFlg;  }

    // Added by Stefan Wiegand at 07.01.2004
    /*! Interface for a buffer that stores performance values indepently from fitness */

    unsigned getNoOfDisplVal() const { return noDisplVal;};
 
    void initDisplVal(unsigned o=0){
      noDisplVal = o;
      displVal.resize(getNoOfDisplVal());
    }

    void setDisplVal     ( unsigned i, double e = 0)  { 
      if(i<displVal.size()) displVal[i] = e;
      else{
	std::cerr << "Individual.h: Individuals error value buffer not correctly initialized!" 
		  << std::endl; exit(0);
      }
    }

    double getDisplVal   ( unsigned i ) const         {
      if(i<displVal.size()) return displVal[i];    
      else{ 
	std::cerr << "Individual.h: Individuals error value buffer not correctly initialized!" 
		  << std::endl; exit(0);
      }
    } 

   //=======================================================================

    bool operator == ( const Individual& ) const;
    bool operator <  ( const Individual& ) const;

    double   getFitness          ( ) const        { return fitness; }

    double   getScaledFitness    ( ) const        { return scaledFitness; }

    void     setScaledFitness    ( double sf )    { scaledFitness = sf; }

    void     setSelProb          ( double sp )    { selProb = sp; }

    double   getSelProb          ( ) const        { return selProb; }

    void     setNumCopies        ( unsigned snc ) { numCopies = snc; }

    unsigned getNumCopies        ( ) const        { return numCopies; }

    void     setEvalFlg          ( bool ef )      { evalFlg = ef; }

    bool     getEvalFlg          ( ) const        { return evalFlg; }

    bool     getFeasible         ( ) const        { return feasible; }

    void     setElitist          ( bool e )       { elitist = e; }

    bool     getElitist          ( ) const        { return elitist; }


    //=======================================================================

    //
    // Added by Marc Toussaint & Stefan Wiegand at 20.11.2002
    //

    /*! Part of PVM-send routine for individuals */
    int pvm_pkind();
   
    /*! Part of PVM-rceive routine for individuals */
    int pvm_upkind();
    

    //=======================================================================

  protected:
    double   fitness;
    double   scaledFitness;
    bool     evalFlg;
    bool     feasible;
    double   selProb;
    unsigned numCopies;
    bool     elitist;
    unsigned age;

    // Added by Stefan Wiegand at 20.11.2002
    unsigned learnTime;

    // Added by Stefan Wiegand at 07.01.2004
    std::vector< double > displVal;
    unsigned noDisplVal;

#ifndef __NO_GENERIC_IOSTREAM
    friend std::ostream& operator << ( std::ostream& os,
				       const Individual& ind ) {
      os << "Individual(" << ind.size( ) << ")\n"
	 << ind.fitness       << '\n'
	 << ind.scaledFitness << '\n'
	 << ind.evalFlg       << '\n'
	 << ind.feasible      << '\n'
	 << ind.selProb       << '\n'
	 << ind.numCopies     << '\n'
	 << ind.elitist       << '\n'
	 << ind.age           << '\n'      
     	 << ind.learnTime     << '\n'    
         << ind.getNoOfDisplVal();
         // added by Stefan Wiegand 07.01.2004
         for(unsigned ii=0; ii<ind.getNoOfDisplVal();ii++)
	   os << '\n' << ind.getDisplVal(ii);
      os << std::endl;
          
      for( unsigned i = 0; i < ind.size( ); ++i )
	  os << '\n' << ind[ i ];
        os << std::endl;
      return os;
    }

    friend std::istream& operator >> ( std::istream& is, Individual& ind ){
      unsigned i, indSize( 0 );
      std::string s, t;
      double e;

      is >> s;
      is.get( );   // skip end of line
      
      if( is.good( ) &&
	  s.substr( 0, 11 ) == "Individual(" &&
	  s.find( ')' ) != std::string::npos ) {

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
 	   >> ind.noDisplVal;
 	// added by Stefan Wiegand 07.01.2004
 	ind.initDisplVal(ind.getNoOfDisplVal());
 	for(i=0; i<ind.getNoOfDisplVal();i++){
 	  is >> e;
 	  ind.setDisplVal(i,e);
 	}
	for( i = 0; i < ind.size( ); ++i )
	  is >> ind[ i ];	
      }
      return is;
    }
#endif // !__NO_GENERIC_IOSTREAM

  friend class Population;
  friend class PopulationMOO;
  friend class IndividualMOO;
};

//===========================================================================

#endif /* !__INDIVIDUAL_H */

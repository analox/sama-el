//===========================================================================
/*!
*
 *  \file Chromosome.h
 *
 *  \author  Martin Kreutz
 *  \date    1995-01-01
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
 *      $RCSfile: Chromosome.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: Chromosome.h,v $
 *      Revision 2.6  2006/04/18 10:23:46  glasmtbl
 *      TG 4/2006: gcc 4.1.0 compatibility ensured
 *
 *      Revision 2.4  2004/05/20 12:23:50  shark-admin
 *      documentation enhanced
 *
 *      Revision 2.3  2004/04/14 12:43:42  shark-admin
 *      new default parameter 'chromswap' added to void crossover(...)-methods that randomly pick crssover-loci of chromosomes
 *
 *      Revision 2.2  2004/04/12 07:43:18  shark-admin
 *      some documentation was changed
 *
 *      Revision 2.1  2004/01/12 11:07:39  shark-admin
 *      documentation modified
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



#ifndef __CHROMOSOME_H
#define __CHROMOSOME_H

#ifdef _WIN32
// disable warning C4786: symbol length > 255 character,
// okay to ignore
#pragma warning(disable: 4786)
// disable warning C4804: suspicious bool comparison
// occurs during instantiation of vector< bool >
#pragma warning(disable: 4804)
#endif

#include <string>
#include <vector>

//
// problems with VC++ 5.0 (can't handle nested scopes)
//
//using namespace std;

//===========================================================================


// forward declarations
class Individual;
class Population;
class IndividualMOO;
class PopulationMOO;


class Chromosome
{
  friend class Individual;
  friend class Population;
  friend class IndividualMOO;
  friend class PopulationMOO;

  public:
    bool sameType( const Chromosome& c ) const
    {
	return strcmp( typeOfAlleles( ), c.typeOfAlleles( ) ) == 0;
    }

    virtual ~Chromosome( ) { }
    virtual const char* typeOfAlleles   ( ) const = 0;
    virtual unsigned    size            ( ) const = 0;
    virtual Chromosome& operator =      ( const Chromosome& ) = 0;

    virtual unsigned    sizeOfAlleles   ( ) const = 0;  
  //
    //=======================================================================
    //
    // string operators
    //
    virtual void        resize          ( unsigned n ) = 0;

    virtual void        duplicate       ( unsigned start,
					  unsigned stop,
					  unsigned dest ) = 0;
    virtual void        invert          ( unsigned start,
					  unsigned stop,
					  unsigned granularity = 1 ) = 0;
    virtual void        transcribe      ( unsigned start,
					  unsigned stop,
					  const Chromosome& chrom ) = 0;
    virtual void        swap            ( unsigned i, unsigned j ) = 0;
    virtual void        shuffle         ( ) = 0;
    virtual void        replace         ( unsigned i,
					  const Chromosome& chrom ) = 0;
    virtual void        insert          ( unsigned i,
					  const Chromosome& chrom ) = 0;
    virtual void        append          ( const Chromosome& chrom ) = 0;
    virtual void        remove          ( unsigned i ) = 0;
    virtual void        remove          ( unsigned from, unsigned to ) = 0;
    virtual void        rotateRight     ( unsigned n = 1 ) = 0;
    virtual void        rotateLeft      ( unsigned n = 1 ) = 0;

#ifndef _WIN32
    //
    // Visual C++ 5.0 can't distinguish between different instantiations of
    // template class vector< T >
    virtual void        crossover
						(
							const Chromosome& dad,
							const Chromosome& mom,
							const std::vector< unsigned >& points
						) = 0;
    virtual void        crossover
						(
							Chromosome& mate,
							const std::vector< unsigned >& points
						) = 0;
#endif // !_WIN32

    virtual void        crossover
						(
							const Chromosome& dad,
							const Chromosome& mom,
							const std::vector< bool >& pos
						) = 0;
    virtual void        crossover
						(
							Chromosome& mate,
							const std::vector< bool >& pos
						) = 0;
    virtual void        crossover
						(
							const Chromosome& dad,
							const Chromosome& mom,
							unsigned npoints,
							unsigned align = 1,
							bool chromswap = 0
						) = 0;
    virtual void        crossover
						(
							Chromosome& mate,
							unsigned npoints,
							unsigned align = 1,
							bool chromswap = 0

						) = 0;
    virtual void        crossover
						(
							const Chromosome& dad,
							const Chromosome& mom,
							const Chromosome& pos
						) = 0;
    virtual void        crossover
						(
							Chromosome& mate,
							const Chromosome& pos
						) = 0;
    virtual void        crossoverUniform
						(
							const Chromosome& dad,
							const Chromosome& mom,
							const std::vector< bool >& pos
						) = 0;
    virtual void        crossoverUniform
						(
							Chromosome& mate,
							const std::vector< bool >& pos
						) = 0;
    virtual void        crossoverUniform
						(
							const Chromosome& dad,
							const Chromosome& mom
						) = 0;
    virtual void        crossoverUniform
						(
							Chromosome& mate
						) = 0;
    virtual void        crossoverUniform
						(
							const Chromosome& dad,
							const Chromosome& mom,
							const Chromosome& pos
						) = 0;
  //  virtual void        crossoverUniform( Chromosome& mate,
  //					  const Chromosome& pos );
    virtual void        recombineDiscrete(const Chromosome& dad,
					  const Chromosome& mom ) = 0;
    virtual void        recombineDiscrete(Chromosome& mate ) = 0;
  //

    virtual bool operator <  ( const Chromosome& c ) const;
    virtual bool operator == ( const Chromosome& c ) const;
 

    // automatically generated
    //virtual bool operator != ( const Chromosome& c ) const;
    //virtual bool operator >  ( const Chromosome& c ) const;
    //virtual bool operator <= ( const Chromosome& c ) const;
    //virtual bool operator >= ( const Chromosome& c ) const;


    //=======================================================================

    //
    // Added by Marc Toussaint & Stefan Wiegand at 20.11.02
    // (interface methods for more externally from Shark defined Chromosomes)

    typedef unsigned uint;
	
    /*! Chromosomes externally defined from the EALib shall'know' how to initialize themselves. */
    virtual void init();

    /*! Chromosomes externally defined from the EALib shall 'know' how to initialize themselves by files. */
    virtual void init(const char* filename);

    /*! Chromosomes externally defined from the EALib shall 'know' how to mutate themselves. */
    virtual void mutate();
 
    /*! Chromosomes externally defined from the EALib shall 'know' the memory loci
        of other Chromosomes located on the same Individual.
	(reasonable, e.g., to read out genes for selfadaptation)
    */
    virtual void registerIndividual(const Individual& i,uint you);

    /*! Shall define a convention for the order of different Chromosomes 
        located on the same Individual and appends all such Chromosomes to the Inividual.
	This method might be called only once per Individual by the 'core' Chromosome that 
	also appends Chromosomes, e.g., of genes dedicated to selfadaptation.
    */
    virtual void appendToIndividual(Individual& i);

    /*! Part of PVM-send routine for all kinds of Chromosomes*/
    virtual int  pvm_pkchrom();

    /*! Part of PVM-receive routine for all kinds of Chromosomes*/
    virtual int  pvm_upkchrom();

    
    //=======================================================================

  protected:
    virtual Chromosome* clone           ( ) const = 0;
    virtual Chromosome* empty           ( ) const = 0;

#ifndef __NO_GENERIC_IOSTREAM
    virtual void        writeTo         ( std::ostream& ) const = 0;
    virtual void        readFrom        ( std::istream& ) = 0;

    friend inline std::istream& operator >> ( std::istream& is, Chromosome& c )
    {
        c.readFrom( is );
		IO_CHECK( is )
		return is;
    }

    friend inline std::ostream& operator << ( std::ostream& os, const Chromosome& c )
    {
        c.writeTo( os );
		IO_CHECK( os )
		return os;
    }
#endif // !__NO_GENERIC_IOSTREAM
};

//===========================================================================

#endif /* !__CHROMOSOME_H */

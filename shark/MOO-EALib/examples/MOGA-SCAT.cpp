/* ======================================================================
 *
 *  File    :  MOGA-SCAT.cpp
 *  Created :  2002-01-23
 *
 *  Description : Sample Program for MOO-ES
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
 *      $RCSfile: MOGA-SCAT.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: MOGA-SCAT.cpp,v $
 *      Revision 2.1  2004/12/14 16:18:32  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.2  2004/12/14 16:03:09  shark-admin
 *      Header corrected
 *
 *      Revision 1.1  2004/11/15 16:21:43  shark-admin
 *      renamed by the addition of the postfix '-SCAT' in order to indicate the dependence on Spread CAT
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


#include "MOO-EALib/PopulationMOO.h"
#include "MOO-EALib/ArchiveMOO.h"
#include "MOO-EALib/TestFunction.h"
#include "Array/Array.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "fstream"
#include "SpreadCAT/SpreadCAT.h"


// main program
int main( void )
{
  // constants
  VarT< unsigned > Seed        ("Seed"       , 1234);
  VarT< bool     > MOOF1       ("MOOF1"      , true);
  VarT< bool     > MOOF2       ("MOOF2"      , false);
  VarT< bool     > MOOF3       ("MOOF3"      , false);
  VarT< bool     > MOOF4       ("MOOF4"      , false);
  VarT< bool     > MOOF5       ("MOOF5"      , false);
  VarT< unsigned > PopSize     ("PopSize"    , 100);
  VarT< unsigned > Dimension   ("Dimension"  , 2);   
  VarT< unsigned > NumOfBits   ("NumOfBits"  , 10);  
  VarT< unsigned > Iterations  ("Iterations" , 1000); 
  VarT< unsigned > DspInterval ("DspInterval", 10);   
  VarT< bool     > Stopping    ("Stopping"   , false);
  VarT< unsigned > NElitists   ("NElitists"  , 0);
  VarT< unsigned > Omega       ("Omega"      , 5);     
  VarT< unsigned > CrossPoints ("CrossPoints", 2);       
  VarT< double   > CrossProb   ("CrossProb"  , 0.6);     
  VarT< double   > FlipProb    ("FlipProb"   , 0.002);
  VarT< bool     > UseGrayCode ("UseGrayCode", true);       
  VarT< double   > Start       ("Start"      , -2);
  VarT< double   > End         ("End"        , 2);
  VarT< unsigned > numArchive  ("numArchive" , 200);
  VarT< double   > Sharing     ("Sharing"    , 0.1);
  VarT< unsigned > FileNumber  ("FileNumber" , 1);
  VarT< bool     > Goldberg    ("Goldberg"   , true);
  VarT< bool     > Fonseca     ("Fonseca"    , false);

  //SpreadCAT
  SpreadCAT mySpread( "MOGA.conf" );
  
  //Other Variable
  const Interval RangeOfValues( Start, End );
  double f1=0, f2=0;
  double PF1[ (unsigned)PopSize ], PF2[ (unsigned)PopSize ];
  double OF1[ (unsigned)PopSize ], OF2[ (unsigned)PopSize ];

  // initialize random number generator
  Rng::seed( Seed );
  // define populations
  PopulationMOO parents   ( PopSize, ChromosomeT< bool >( Dimension * NumOfBits ) );
  PopulationMOO offsprings( PopSize, ChromosomeT< bool >( Dimension * NumOfBits ) );
  // Set Minimization Task
  parents   .setMinimize( );
  offsprings.setMinimize( );
  // Set No of Objective funstions
  parents   .setNoOfObj( 2 );
  offsprings.setNoOfObj( 2 );
  // Archive
  ArchiveMOO archive( numArchive );
  archive.minimize( );
  // Temporary chromosomes for decoding
  ChromosomeT< double > dblchrom;
  // initialize all chromosomes of parent population
  for ( unsigned i = 0; i < parents.size( ); ++i ){
    dynamic_cast< ChromosomeT< bool >& >( parents[ i ][ 0 ] ).initialize( );
  }
  // evaluate parents (only needed for elitist strategy)
  if ( NElitists > 0 ){
    for ( unsigned i = 0; i < parents.size( ); ++i ) {
      dblchrom.decodeBinary( parents[ i ][ 0 ], RangeOfValues, NumOfBits, UseGrayCode );
      if ( MOOF1 ){
        f1 = SphereF1( dblchrom );
        f2 = SphereF2( dblchrom );
      }
      else if ( MOOF2 ){
	f1 = DebConvexF1( dblchrom );
	f2 = DebConvexF2( dblchrom );
      }
      else if ( MOOF3 ){
        f1 = DebConcaveF1( dblchrom );
	f2 = DebConcaveF2( dblchrom );
      }
      else if ( MOOF4 ){
	f1 = DebDiscreteF1( dblchrom );
	f2 = DebDiscreteF2( dblchrom );
      }
      else if ( MOOF5 ){
	f1 = FonsecaConcaveF1( dblchrom );
	f2 = FonsecaConcaveF2( dblchrom );
      }
      parents[ i ].setMOOFitnessValues( f1, f2 );
      PF1[ i ] = f1;
      PF2[ i ] = f2;
    }
    parents.NicheCountPFN2( Sharing );
    if ( Goldberg ){
      parents.MOGAGoldbergRank( );
    }
    else if ( Fonseca ){
      parents.MOGAFonsecaRank( );
    }
    parents.MOORankToFitness( );
    parents.SelectProbMichalewicz( );
    parents.SharingSelProb( );
  }
  // iterate
  for ( unsigned t = 1; t < Iterations+1; ++t ) {
    std::cout << "Generations : " << t << std::endl;
    // copy parents to offsprings
    offsprings = parents;
    // recombine by crossing over two parents
    for ( unsigned i = 0; i < offsprings.size( )-1; i += 2 ){
      if ( Rng::coinToss( CrossProb ) ){
        offsprings[ i ][ 0 ].crossover( offsprings[ i+1 ][ 0 ], CrossPoints );
      }
    }
    // mutate by flipping bits
    for ( unsigned i = 0; i < offsprings.size( ); ++i ){
      dynamic_cast< ChromosomeT< bool >& >( offsprings[ i ][ 0 ] ).flip( FlipProb );
    }
    // evaluate objective function
    for ( unsigned i = 0; i < offsprings.size( ); ++i ) {
      dblchrom.decodeBinary( offsprings[ i ][ 0 ], RangeOfValues, NumOfBits, UseGrayCode );
      if ( MOOF1 ){
        f1 = SphereF1( dblchrom );
        f2 = SphereF2( dblchrom );
      }
      else if ( MOOF2 ){
	f1 = DebConvexF1( dblchrom );
	f2 = DebConvexF2( dblchrom );
      }
      else if ( MOOF3 ){
        f1 = DebConcaveF1( dblchrom );
	f2 = DebConcaveF2( dblchrom );
      }
      else if ( MOOF4 ){
	f1 = DebDiscreteF1( dblchrom );
	f2 = DebDiscreteF2( dblchrom );
      }
      else if ( MOOF5 ){
	f1 = FonsecaConcaveF1( dblchrom );
	f2 = FonsecaConcaveF2( dblchrom );
      }
      offsprings[ i ].setMOOFitnessValues( f1, f2 );
      OF1[ i ] = f1;
      OF2[ i ] = f2;
    }
    offsprings.NicheCountPFN2( Sharing );
    if ( Goldberg ){
      offsprings.MOGAGoldbergRank( );
    }
    else if ( Fonseca ){
      offsprings.MOGAFonsecaRank( );
    }
    offsprings.MOORankToFitness( );
    offsprings.SelectProbMichalewicz( );
    offsprings.SharingSelProb( );
    // Archive
    for ( unsigned i = 0; i < offsprings.size( ); i++ ){
      int dominateIA = archive.Dominate( offsprings[ i ] );
      if ( dominateIA >= 4 ){
	archive.cleanArchive( );
	archive.addArchive( offsprings[ i ] );
      }
      else if ( dominateIA == 3 ){
	archive.delDominateArchive( offsprings[ i ] );
	archive.addArchive( offsprings[ i ] );
      }
      else if ( dominateIA == 2 ){
	if ( archive.getCapacity( ) > 0 ){
	  archive.addArchive( offsprings[ i ] );
	}
	else {
	  double mindisIA = archive.distanceOnFitness( offsprings[ i ] );
	  double mindisAA = archive.minDistanceOnFitness( );
          if ( mindisIA > mindisAA ){
	    archive.delSharingWorst( );
	    archive.addArchive( offsprings[ i ] );
	  }
	}
      }
    }

    // selection
    parents.selectElitists( offsprings, NElitists );
    parents.SelectByRoulette( offsprings, NElitists );
    // Evaluate again ( Selection Copy does not work well )
    for ( unsigned i = 0; i < parents.size( ); ++i ) {
      dblchrom.decodeBinary( parents[ i ][ 0 ], RangeOfValues, NumOfBits, UseGrayCode );
      if ( MOOF1 ){
        f1 = SphereF1( dblchrom );
        f2 = SphereF2( dblchrom );
      }
      else if ( MOOF2 ){
	f1 = DebConvexF1( dblchrom );
	f2 = DebConvexF2( dblchrom );
      }
      else if ( MOOF3 ){
        f1 = DebConcaveF1( dblchrom );
	f2 = DebConcaveF2( dblchrom );
      }
      else if ( MOOF4 ){
	f1 = DebDiscreteF1( dblchrom );
	f2 = DebDiscreteF2( dblchrom );
      }
      else if ( MOOF5 ){
	f1 = FonsecaConcaveF1( dblchrom );
	f2 = FonsecaConcaveF2( dblchrom );
      }
      parents[ i ].setMOOFitnessValues( f1, f2 );
      PF1[ i ] = f1;
      PF2[ i ] = f2;
    }
    parents.NicheCountPFN2( Sharing );
    if ( Goldberg ){
      parents.MOGAGoldbergRank( );
    }
    else if ( Fonseca ){
      parents.MOGAFonsecaRank( );
    }
    parents.MOORankToFitness( );
    parents.SelectProbMichalewicz( );
    parents.SharingSelProb( );

    if ( t % DspInterval == 0 ){
      std::cout << "Archive Size = " << archive.size( ) << std::endl;
      if ( Stopping ){
	std::cout << "\nHit return key" << std::endl;
	fgetc( stdin );
      }
    }
  }
    
  char filename[ 80 ];
  unsigned filenumber = FileNumber;
  sprintf( filename, "Archive%i.txt", filenumber );
  archive.saveArchive( filename );
  std::cout << "\nThe result is stored in " << filename << "\n" << std::endl;
  std::cout << "Number of solutions" << std::endl;
  std::cout << "Number of objectives" << std::endl;
  std::cout << "Solution 1 (f1,f2,...)" << std::endl;
  std::cout << "Solution 2 (f1,f2,...)\n" << std::endl;
  
  std::cout << "\nHit return key" << std::endl;
  
  fgetc(stdin);

  return 0;

}

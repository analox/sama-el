/* ======================================================================
 *
 *  File    :  DWA-SCAT.cpp
 *  Created :  2001-12-16
 *
 *  Description : Sample Program for MOO-ES ( Standard )
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
 *      $RCSfile: DWA-SCAT.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: DWA-SCAT.cpp,v $
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


#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>
#include <MOO-EALib/TestFunction.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Array/Array.h>
#include <SpreadCAT/SpreadCAT.h>

//=======================================================================
// main program
//=======================================================================
int main( void )
{

  //Constants & Variables
  VarT< bool     > Standard1    ("Standard1",TRUE);
  VarT< bool     > Standard2    ("Standard2",FALSE);
  VarT< bool     > Derandomize  ("Derandomize",FALSE);
  VarT< bool     > Rotate       ("Rotate",FALSE);
  VarT< bool     > GSA          ("GSA",FALSE);
  VarT< bool     > IDA          ("IDA",FALSE);
  VarT< bool     > CMA          ("CMA",FALSE);
  VarT< unsigned > baseSize     ("baseSize",1350);
  VarT< unsigned > Seed         ("Seed", 1234);
  VarT< unsigned > MU           ("MU",15);
  VarT< unsigned > LAMBDA       ("LAMBDA",100);
  VarT< bool     > PLUSSTRATEGY ("PLUSSTRATEGY",FALSE);
  VarT< bool     > MOOF1        ("MOOF1",TRUE);
  VarT< bool     > MOOF2        ("MOOF2",FALSE);
  VarT< bool     > MOOF3        ("MOOF3",FALSE);
  VarT< bool     > MOOF4        ("MOOF4",FALSE);
  VarT< bool     > MOOF5        ("MOOF5",FALSE);
  VarT< bool     > USEDWA       ("USEDWA",TRUE);
  VarT< bool     > USEBWA       ("USEBWA",FALSE);
  VarT< unsigned > FREQUENCY    ("FREQUENCY",500);
  VarT< double   > PHASE        ("PHASE",0);
  VarT< bool     > Recombine    ("Recombine",FALSE); 
  VarT< unsigned > Dimension    ("Dimension",2);
  VarT< unsigned > Iterations   ("Iterations",2000);
  VarT< bool     > CheckSigma   ("CheckSigma",FALSE); 
  VarT< double   > SigmaLower   ("SigmaLower",0.004);
  VarT< unsigned > numArchive   ("numArchive",200);
  VarT< double   > NICHE        ("NICHE",0.05);
  VarT< double   > MinInit      ("MinInit",-2);  
  VarT< double   > MaxInit      ("MaxInit",2); 
  VarT< double   > SigmaMin     ("SigmaMin",1);
  VarT< double   > SigmaMax     ("SigmaMax",1);
  VarT< unsigned > Interval     ("Interval",100);
  VarT< bool     > Stopping     ("Stopping",FALSE);
  VarT< unsigned > FileName     ("FileName",1);

  unsigned NSigma;
  unsigned NELITE;
  std::vector< double > weight( 2 );
  int      domination;
  double   mindistanceIA, mindistanceAA;

  //char MatlabCommand[200];
  //double Var1, Var2, Var3, Var4;

  //SpreadCAT
  SpreadCAT mySpread( "DWA.conf" );
  //Other Variables
  if ( Standard1 ){
    NSigma = Dimension;
  }
  else if ( Standard2 ){
    NSigma = Dimension;
  }
  else if ( Derandomize ){
    std::cout << "Not defined" << std::endl;
    exit(-1);
  }
  else if ( Rotate ){
    NSigma = (unsigned)( 0.5 * ( 1 + Dimension ) * Dimension );
  }
  else if ( GSA ){
    NSigma = 1 + baseSize * Dimension;
  }
  else if ( IDA ){
    NSigma = 2 + 3 * Dimension;
  }
  else if ( CMA ){
    NSigma = 1 + ( 2 + Dimension ) * Dimension;
  }
  else {
    exit(-1);
  }
  double           XP[ (unsigned)MU ], YP[ (unsigned)MU ];
  double           XO[ (unsigned)LAMBDA ], YO[ (unsigned)LAMBDA ];
  double           XA[ (unsigned)numArchive ], YA[ (unsigned)numArchive ];
  double           f1=0, f2=0;
  // First Weight
  if ( USEDWA ){
    weight[ 0 ] = 0.5 + 0.5 * sin( PHASE * 3.141592 );
    weight[ 1 ] = 1.0 - weight[ 0 ];
  }
  else if ( USEBWA ){
    if ( (unsigned)PHASE % 2 == 0 ){
      weight[ 0 ] = 1.0;
    }
    else {
      weight[ 0 ] = 0.0;
    }
    weight[ 1 ] = 1.0 - weight[ 0 ];
  }
  //Random Seed
  Rng::seed( Seed );
  //Define Populations
  PopulationMOO parents   ( MU,     
                ChromosomeT< double >( Dimension ), ChromosomeT< double >( NSigma ));
  PopulationMOO offsprings( LAMBDA, 
                ChromosomeT< double >( Dimension ), ChromosomeT< double >( NSigma ));
  //Set the Number of Objective Functions
  parents   .setNoOfObj( 2 );
  offsprings.setNoOfObj( 2 );
  //Task ( Minimize or Maximize )
  parents.   setMinimize( );
  offsprings.setMinimize( );
  //Archive
  ArchiveMOO archive( (unsigned)numArchive );
  archive.minimize( );
  //Initialize Parent Population
  for ( unsigned i = 0; i < MU; i++ ){
    dynamic_cast< ChromosomeT< double >& >( parents[i][0] ).initialize( MinInit, MaxInit );
    if ( Standard1 ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initialize( SigmaMin, SigmaMax );
    }
    else if ( Standard2 ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initialize( SigmaMin, SigmaMax );
    }
    else if ( Derandomize ){
      std::cout << "Not defined" << std::endl;
      exit(-1);
    }
    else if ( Rotate ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeRotate( SigmaMin, SigmaMax );
    }
    else if ( GSA ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeGSA( SigmaMin, SigmaMax, baseSize );
    }
    else if ( IDA ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeIDA( SigmaMin, SigmaMax );
    }
    else if ( CMA ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeCMA( SigmaMin, SigmaMax );
    }
    else {
      exit(-1);
    }
  }
  // PLUS Strategy
  if ( PLUSSTRATEGY ){
    NELITE = MU;
  }
  else {
    NELITE = 0;
  }
  //Evaluate Parents ( only needed for elitist strategy )
  if ( PLUSSTRATEGY ){
    for ( unsigned i = 0; i < MU; i++ ){
      if ( MOOF1 ){
        f1 = SphereF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = SphereF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF2 ){
        f1 = DebConvexF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = DebConvexF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF3 ){
        f1 = DebConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = DebConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF4 ){
        f1 = DebDiscreteF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = DebDiscreteF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF5 ){
        f1 = FonsecaConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = FonsecaConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      parents[i].setMOOFitnessValues( f1, f2 );
      XP[ i ] = f1;
      YP[ i ] = f2;
    }
    parents.aggregation( weight );
  }
  //Iterate
  for( unsigned t = 1; t < Iterations + 1; t++ ) {
    //Weight
    if ( USEDWA ){
      weight[ 0 ] = 0.5 + 0.5 * sin( ( PHASE + (double)(t-1)/(double)FREQUENCY ) * 3.141592 );
      weight[ 1 ] = 1.0 - weight[ 0 ];
    }
    else if ( USEBWA ){
      if ( (unsigned)(PHASE+(double)(t-1)/(double)FREQUENCY) % 2 == 0 ){
        weight[ 0 ] = 1.0;
      }
      else {
        weight[ 0 ] = 0.0;
      }
      weight[ 1 ] = 1.0 - weight[ 0 ];
    }
    //Print Generation
    std::cout << "Generation  = " << t;
    std::cout << "   archive  = " << archive.size( );
    std::cout << "   w1 = " << weight[ 0 ];
    std::cout << "   w2 = " << weight[ 1 ] << std::endl;
    // Calculate Fitness in Parents
    parents.aggregation( weight );
    //Generate New Offsprings
    for ( unsigned i = 0; i < LAMBDA; i++ ){
      IndividualMOO& pap = parents.random( );
      IndividualMOO& mom = parents.random( );
      ChromosomeT< double >& objvar = dynamic_cast< ChromosomeT< double >& >( offsprings[i][0] );
      ChromosomeT< double >& sigma  = dynamic_cast< ChromosomeT< double >& >( offsprings[i][1] );
      //Recombination & Mutation
      if ( Standard1 ){
        if ( Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
          sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
        else {
	  for ( unsigned j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( unsigned j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
        double tau0 = 1.0 / sqrt( 2.0 * (double)Dimension );
        double tau1 = 1.0 / sqrt( 2.0 * sqrt( (double)Dimension ));
        sigma.mutateLogNormal( tau0, tau1 );
        if (CheckSigma){
          for ( unsigned j = 0; j < sigma.size( ); j++ ){
            if ( sigma[ j ] < SigmaLower * fabs( objvar[ j ] ) ){
              sigma[ j ] = SigmaLower * fabs( objvar[ j ] );
            }
          }
        }
        objvar.mutateNormal( sigma, true );
      }
      else if ( Standard2 ){
        if ( Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
          sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
        else {
	  for ( unsigned j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( unsigned j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
        double xi_prob = 0.5;
        sigma.mutateMSR( xi_prob );
        for ( unsigned j = 0; j < sigma.size( ); j++ ){
	  sigma[ j ] = sigma[ j ] / sqrt( sigma.size( ) );
        }
        if (CheckSigma){
          for ( unsigned j = 0; j < sigma.size( ); j++ ){
            if ( sigma[ j ] < SigmaLower * fabs( objvar[ j ] ) ){
              sigma[ j ] = SigmaLower * fabs( objvar[ j ] );
            }
          }
        }
        objvar.mutateNormal( sigma, true );
        for ( unsigned j = 0; j < sigma.size( ); j++ ){
	  sigma[ j ] = sigma[ j ] * sqrt( sigma.size( ) );
        }
      }
      else if ( Derandomize ){
        std::cout << "Not defined" << std::endl;
        exit(-1);
      }
      else if ( Rotate ){
        if ( Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
          sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
        else {
	  for ( unsigned j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( unsigned j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
        objvar.mutateRotate( sigma );
      }
      else if ( GSA ){
        if ( Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
          sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
        else {
	  for ( unsigned j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( unsigned j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
        objvar.mutateGSA( sigma );
      }
      else if ( IDA ){
        if ( Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
          sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
        else {
	  for ( unsigned j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( unsigned j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
        objvar.mutateIDA( sigma );
      }
      else if ( CMA ){
        if ( Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
          sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
        else {
	  for ( unsigned j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
          }
          for ( unsigned j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
        objvar.mutateCMA( sigma );
      }
      else {
        exit(-1);
      }    
    }
    //Evaluate Objective Function
    for ( unsigned j = 0; j < LAMBDA; j++ ){
      if ( MOOF1 ){
        f1 = SphereF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
        f2 = SphereF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
      }
      if ( MOOF2 ){
        f1 = DebConvexF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
        f2 = DebConvexF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
      }
      if ( MOOF3 ){
        f1 = DebConcaveF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
        f2 = DebConcaveF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
      }
      if ( MOOF4 ){
        f1 = DebDiscreteF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
        f2 = DebDiscreteF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
      }
      if ( MOOF5 ){
        f1 = FonsecaConcaveF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
        f2 = FonsecaConcaveF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
      }
      offsprings[j].setMOOFitnessValues( f1, f2 );
      XO[ j ] = f1;
      YO[ j ] = f2;
    }
    offsprings.aggregation( weight );
   
    // Selection ( Mu, Labmda ) or ( Mu + Lambda )
    parents.selectMuLambda( offsprings, NELITE );


    for ( unsigned i = 0; i < MU; i++ ){
      if ( MOOF1 ){
        f1 = SphereF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = SphereF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF2 ){
        f1 = DebConvexF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = DebConvexF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF3 ){
        f1 = DebConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = DebConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF4 ){
        f1 = DebDiscreteF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = DebDiscreteF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      if ( MOOF5 ){
        f1 = FonsecaConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
        f2 = FonsecaConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
      }
      parents[i].setMOOFitnessValues( f1, f2 );
      XP[ i ] = f1;
      YP[ i ] = f2;
    }
    parents.aggregation( weight );

    // Store Parent Data
    for ( unsigned j = 0; j < parents.size( ); j++ ){
      XP[ j ] = parents[ j ].getMOOFitness( 0 );
      YP[ j ] = parents[ j ].getMOOFitness( 1 );
    }

    // Archive
    for ( unsigned j = 0; j < offsprings.size( ); j++ ){
      domination = archive.Dominate( offsprings[ j ] );
      mindistanceIA = archive.distanceOnFitness( offsprings[ j ] );
      if ( domination > 2 ){
        archive.delDominateArchive( offsprings[ j ] );
	archive.addArchive( offsprings[ j ] );
      }
      else if ( domination == 2 ){
	if ( mindistanceIA > NICHE ){
          if ( archive.getCapacity( ) > 0 ){
	    archive.addArchive( offsprings[ j ] );
	  }
	  else {
            mindistanceAA = archive.minDistanceOnFitness( );
            if ( mindistanceIA > mindistanceAA ){
              archive.delSharingWorst( );
	      archive.addArchive( offsprings[ j ] );
	    }
	    else {
	    }
	  }
	}
	else {
	}
      }
      else {
      }
    }

    // Store Archive Data
    for ( unsigned j = 0; j < archive.size( ); j++ ){
      XA[ j ] = archive.readArchive( j ).getMOOFitness( 0 );
      YA[ j ] = archive.readArchive( j ).getMOOFitness( 1 );
    }

    if ( t % Interval == 0 ){
      if ( Stopping ){
        std::cout << "\n Hit Return Key for Continue" << std::endl;
        fgetc(stdin);
      }
    }
  }


  char dummy[ 80 ];
  unsigned filenumber = FileName;
  sprintf( dummy, "Archive%i.txt", filenumber );
  archive.saveArchive( dummy );
  std::cout << "\nThe result is stored in " << dummy << "\n" << std::endl;
  std::cout << "Number of Solutions" << std::endl;
  std::cout << "Number of Objectives" << std::endl;
  std::cout << "Solution 1 (f1,f2,...)" << std::endl;
  std::cout << "Solution 2 (f1,f2,...)\n" << std::endl;

  std::cout << "Hit Return Key for End" << std::endl;
  fgetc(stdin);

  return 0;

}


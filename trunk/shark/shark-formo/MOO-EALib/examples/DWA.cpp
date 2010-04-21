/* ======================================================================
 *
 *  File    :  DWA.cpp
 *  Created :  2004-11-15
 *
 *  Description : Sample Program for MOO-ES ( Standard )
 *  Author  : Stefan Roth 
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
 *      $RCSfile: DWA.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: DWA.cpp,v $
 *      Revision 2.9  2005/01/14 15:04:38  shark-admin
 *      uint to double cast
 *
 *      Revision 2.8  2004/12/30 14:12:41  shark-admin
 *      comment added to MyParam::readParam()
 *
 *      Revision 2.7  2004/12/15 15:42:46  shark-admin
 *      pragma warning disabled
 *
 *      Revision 2.6  2004/12/15 13:56:28  shark-admin
 *      switch for fitness functions
 *
 *      Revision 2.5  2004/12/14 17:13:25  shark-admin
 *      header and body formated
 *
 *      Revision 2.4  2004/11/15 16:23:59  shark-admin
 *      examples now independent of SpreadCAT - compatibility with vc++ ensured
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
/* This file is derived from DWA-SCAT.cpp which was written by  
   Tatsuya Okabe <tatsuya.okabe@honda-ri.de>.
**********************************************************/

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

#include <iostream> 
#include <iomanip>
#include <string>
#include <fstream>

#include <MOO-EALib/PopulationMOO.h>
#include <MOO-EALib/ArchiveMOO.h>
#include <MOO-EALib/TestFunction.h>
#include <Array/Array.h>

#include <ReClaM/Params.h>
//#include <SpreadCAT/SpreadCAT.h>

// My own derived class for managing my configuration files:
//
class MyParams : public Params 
{
public:
	
  MyParams( int argc, char **argv ) : Params( argc, argv ) { }
	
  ~MyParams( ) { }
	
  void readParams()
  {
    // Call the main program with parameters "-conf [filename]"
    if ( scanFrom(confFile.c_str()) ) 
      {	
	std::cout << "Name of the configuration file: " << confFile << std::endl;
      } 
    else
      {
	std::cerr << "No valid configuration file given! (" << confFile  << ")" << std::endl;
	exit(0);
      }
  }
	
  void io( std::istream& is, std::ostream& os, FileUtil::iotype type )
  {
    FileUtil::io( is, os, "Standard1"    , Standard1, true, type);
    FileUtil::io( is, os, "Standard2"    , Standard2, false, type);
    FileUtil::io( is, os, "Derandomize"  , Derandomize, false, type);
    FileUtil::io( is, os, "Rotate"       , Rotate, false, type);
    FileUtil::io( is, os, "GSA"          , GSA, false, type);
    FileUtil::io( is, os, "IDA"          , IDA, false, type);
    FileUtil::io( is, os, "CMA"          , CMA, false, type);
    FileUtil::io( is, os, "baseSize"     , baseSize, 1350u, type);
    FileUtil::io( is, os, "Seed"         , Seed,  1234u, type);
    FileUtil::io( is, os, "MU"           , MU, 15u, type);
    FileUtil::io( is, os, "LAMBDA"       , LAMBDA, 100u, type);
    FileUtil::io( is, os, "PLUSSTRATEGY" , PLUSSTRATEGY, false, type);
    FileUtil::io( is, os, "MOOF"         , MOOF, 1u, type);
    FileUtil::io( is, os, "USEDWA"       , USEDWA, true, type);
    FileUtil::io( is, os, "USEBWA"       , USEBWA, false, type);
    FileUtil::io( is, os, "FREQUENCY"    , FREQUENCY, 500u, type);
    FileUtil::io( is, os, "PHASE"        , PHASE, .0, type);
    FileUtil::io( is, os, "Recombine"    , Recombine, false, type); 
    FileUtil::io( is, os, "Dimension"    , Dimension, 2u, type);
    FileUtil::io( is, os, "Iterations"   , Iterations, 2000u, type);
    FileUtil::io( is, os, "CheckSigma"   , CheckSigma, false, type); 
    FileUtil::io( is, os, "SigmaLower"   , SigmaLower, 0.004, type);
    FileUtil::io( is, os, "numArchive"   , numArchive, 200u, type);
    FileUtil::io( is, os, "NICHE"        , NICHE, 0.05, type);
    FileUtil::io( is, os, "MinInit"      , MinInit, -2., type);  
    FileUtil::io( is, os, "MaxInit"      , MaxInit, 2., type); 
    FileUtil::io( is, os, "SigmaMin"     , SigmaMin, 1., type);
    FileUtil::io( is, os, "SigmaMax"     , SigmaMax, 1., type);
    FileUtil::io( is, os, "Interval"     , Interval, 100u, type);
    FileUtil::io( is, os, "Stopping"     , Stopping, false, type);
    FileUtil::io( is, os, "FileName"     , FileName, 1u, type);
  }
	
  // Method to show the current content of the class variable: 
  void monitor( )
  {
    std::cout << "Standard1 "	<< Standard1	<< std::endl;
    std::cout << "Standard2 "	<< Standard2	<< std::endl;
    std::cout << "Derandomize "	<< Derandomize	<< std::endl;
    std::cout << "Rotate "	<< Rotate	<< std::endl;
    std::cout << "GSA "		<< GSA		<< std::endl;
    std::cout << "IDA "		<< IDA		<< std::endl;
    std::cout << "CMA "		<< CMA		<< std::endl;
    std::cout << "baseSize "	<< baseSize	<< std::endl;
    std::cout << "Seed  "	<< Seed		<< std::endl;
    std::cout << "MU "		<< MU		<< std::endl;
    std::cout << "LAMBDA "	<< LAMBDA	<< std::endl;
    std::cout << "PLUSSTRATEGY "<< PLUSSTRATEGY	<< std::endl;
    std::cout << "MOOF "	<< MOOF 	<< std::endl;
    std::cout << "USEDWA "	<< USEDWA	<< std::endl;
    std::cout << "USEBWA "	<< USEBWA	<< std::endl;
    std::cout << "FREQUENCY "	<< FREQUENCY	<< std::endl;
    std::cout << "PHASE "	<< PHASE	<< std::endl;
    std::cout << "Recombine "	<< Recombine	<< std::endl;
    std::cout << "Dimension "	<< Dimension	<< std::endl;
    std::cout << "Iterations "	<< Iterations	<< std::endl;
    std::cout << "CheckSigma "	<< CheckSigma	<< std::endl;
    std::cout << "SigmaLower "	<< SigmaLower	<< std::endl;
    std::cout << "numArchive "	<< numArchive	<< std::endl;
    std::cout << "NICHE "	<< NICHE	<< std::endl;
    std::cout << "MinInit "	<< MinInit	<< std::endl;
    std::cout << "MaxInit "	<< MaxInit	<< std::endl;
    std::cout << "SigmaMin "	<< SigmaMin	<< std::endl;
    std::cout << "SigmaMax "	<< SigmaMax	<< std::endl;
    std::cout << "Interval "	<< Interval	<< std::endl;
    std::cout << "Stopping "	<< Stopping	<< std::endl;
    std::cout << "FileName "	<< FileName	<< std::endl;
  }
	
  bool Standard1;
  bool Standard2;
  bool Derandomize;
  bool Rotate;
  bool GSA;
  bool IDA;
  bool CMA;
  unsigned baseSize;
  unsigned Seed;
  unsigned MU;
  unsigned LAMBDA;
  bool PLUSSTRATEGY;
  unsigned MOOF;
  bool USEDWA;
  bool USEBWA;
  unsigned FREQUENCY;
  double PHASE;
  bool Recombine;
  unsigned Dimension;
  unsigned Iterations;
  bool CheckSigma;
  double SigmaLower;
  unsigned numArchive;
  double NICHE;
  double MinInit;
  double MaxInit;
  double SigmaMin;
  double SigmaMax;
  unsigned Interval;
  bool Stopping;
  unsigned FileName ;
};

// main program
int main (int argc, char* argv[])
{
  unsigned i,j;

  MyParams param( argc, argv );
   	
  param.setDefault( );
  //param.readParams();
  param.monitor();

  unsigned NSigma;
  unsigned NELITE;
  std::vector< double > weight( 2 );
  int      domination;
  double   mindistanceIA, mindistanceAA;
 
  //char MatlabCommand[200];
  //double Var1, Var2, Var3, Var4;

  //Other Variables
  if ( param.Standard1 ){
    NSigma = param.Dimension;
  }
  else if ( param.Standard2 ){
    NSigma = param.Dimension;
  }
  else if ( param.Derandomize ){
    std::cout << "Not defined" << std::endl;
    exit(-1);
  }
  else if ( param.Rotate ){
    NSigma = (unsigned)( 0.5 * ( 1 + param.Dimension ) * param.Dimension );
  }
  else if ( param.GSA ){
    NSigma = 1 + param.baseSize * param.Dimension;
  }
  else if ( param.IDA ){
    NSigma = 2 + 3 * param.Dimension;
  }
  else if ( param.CMA ){
    NSigma = 1 + ( 2 + param.Dimension ) * param.Dimension;
  }
  else {
    exit(-1);
  }

  double *XP = new double[ (unsigned)param.MU ];
  double *YP = new double[ (unsigned)param.MU ];	
  double *XO = new double[ (unsigned)param.LAMBDA ];
  double *YO = new double[ (unsigned)param.LAMBDA ];	
  double *XA = new double[ (unsigned)param.numArchive ];
  double *YA = new double[ (unsigned)param.numArchive ];	
	
  double           f1=0, f2=0;
  // First Weight
  if ( param.USEDWA ){
    weight[ 0 ] = 0.5 + 0.5 * sin( param.PHASE * 3.141592 );
    weight[ 1 ] = 1.0 - weight[ 0 ];
  }
  else if ( param.USEBWA ){
    if ( (unsigned)param.PHASE % 2 == 0 ){
      weight[ 0 ] = 1.0;
    }
    else {
      weight[ 0 ] = 0.0;
    }
    weight[ 1 ] = 1.0 - weight[ 0 ];
  }
  //Random Seed
  Rng::seed( param.Seed );

  //Define Populations
  PopulationMOO parents   ( param.MU,     
			    ChromosomeT< double >( param.Dimension ), ChromosomeT< double >( NSigma ));
  PopulationMOO offsprings( param.LAMBDA, 
			    ChromosomeT< double >( param.Dimension ), ChromosomeT< double >( NSigma ));

  //Set the Number of Objective Functions
  parents   .setNoOfObj( 2 );
  offsprings.setNoOfObj( 2 );

  //Task ( Minimize or Maximize )
  parents.   setMinimize( );
  offsprings.setMinimize( );

  //Archive
  ArchiveMOO archive( (unsigned)param.numArchive );
  archive.minimize( );

  //Initialize Parent Population
  for ( i = 0; i < param.MU; i++ ){
    dynamic_cast< ChromosomeT< double >& >( parents[i][0] ).initialize( param.MinInit, param.MaxInit );
    if ( param.Standard1 ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initialize( param.SigmaMin, param.SigmaMax );
    }
    else if ( param.Standard2 ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initialize( param.SigmaMin, param.SigmaMax );
    }
    else if ( param.Derandomize ){
      std::cout << "Not defined" << std::endl;
      exit(-1);
    }
    else if ( param.Rotate ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeRotate( param.SigmaMin, param.SigmaMax );
    }
    else if ( param.GSA ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeGSA( param.SigmaMin, param.SigmaMax, param.baseSize );
    }
    else if ( param.IDA ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeIDA( param.SigmaMin, param.SigmaMax );
    }
    else if ( param.CMA ){
      dynamic_cast< ChromosomeT< double >& >( parents[i][1] ).initializeCMA( param.SigmaMin, param.SigmaMax );
    }
    else {
      exit(-1);
    }
  }

  // PLUS Strategy
  if ( param.PLUSSTRATEGY ){
    NELITE = param.MU;
  }
  else {
    NELITE = 0;
  }

  //Evaluate Parents ( only needed for elitist strategy )
  if ( param.PLUSSTRATEGY ){
    for ( i = 0; i < param.MU; i++ ){
      switch(param.MOOF)
	{
	case 1:
	  f1 = SphereF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = SphereF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 2:
	  f1 = DebConvexF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = DebConvexF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 3:
	  f1 = DebConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = DebConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 4:
	  f1 = DebDiscreteF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = DebDiscreteF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 5:
	  f1 = FonsecaConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = FonsecaConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	default:
	  std::cerr << "Specify fitness function. Current setting of MOOF (1-5):" << param.MOOF << std::endl;
	  break;
	}
      parents[i].setMOOFitnessValues( f1, f2 );
      XP[ i ] = f1;
      YP[ i ] = f2;
    }
    parents.aggregation( weight );
  }

  //Iterate
  for( unsigned t = 1; t < param.Iterations + 1; t++ ) {
    //Weight
    if ( param.USEDWA ){
      weight[ 0 ] = 0.5 + 0.5 * sin( ( param.PHASE + (double)(t-1)/(double)param.FREQUENCY ) * 3.141592 );
      weight[ 1 ] = 1.0 - weight[ 0 ];
    }
    else if ( param.USEBWA ){
      if ( (unsigned)(param.PHASE+(double)(t-1)/(double)param.FREQUENCY) % 2 == 0 ){
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
    for ( i = 0; i < param.LAMBDA; i++ ){
      IndividualMOO& pap = parents.random( );
      IndividualMOO& mom = parents.random( );
      ChromosomeT< double >& objvar = dynamic_cast< ChromosomeT< double >& >( offsprings[i][0] );
      ChromosomeT< double >& sigma  = dynamic_cast< ChromosomeT< double >& >( offsprings[i][1] );

      //Recombination & Mutation
      if ( param.Standard1 ){
	if ( param.Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
	  sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
	else {
	  for ( j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
	double tau0 = 1.0 / sqrt( 2.0 * (double)param.Dimension );
	double tau1 = 1.0 / sqrt( 2.0 * sqrt( (double)param.Dimension ));
	sigma.mutateLogNormal( tau0, tau1 );
	if (param.CheckSigma){
	  for ( j = 0; j < sigma.size( ); j++ ){
	    if ( sigma[ j ] < param.SigmaLower * fabs( objvar[ j ] ) ){
	      sigma[ j ] = param.SigmaLower * fabs( objvar[ j ] );
	    }
	  }
	}
	objvar.mutateNormal( sigma, true );
      }
      else if ( param.Standard2 ){
	if ( param.Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
	  sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
	else {
	  for ( j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
	double xi_prob = 0.5;
	sigma.mutateMSR( xi_prob );
	for ( j = 0; j < sigma.size( ); j++ ){
	  sigma[ j ] = sigma[ j ] / sqrt((double)sigma.size( ) );
	}
	if (param.CheckSigma){
	  for ( j = 0; j < sigma.size( ); j++ ){
	    if ( sigma[ j ] < param.SigmaLower * fabs( objvar[ j ] ) ){
	      sigma[ j ] = param.SigmaLower * fabs( objvar[ j ] );
	    }
	  }
	}
	objvar.mutateNormal( sigma, true );
	for ( j = 0; j < sigma.size( ); j++ ){
	  sigma[ j ] = sigma[ j ] * sqrt( sigma.size( ) );
	}
      }
      else if ( param.Derandomize ){
	std::cout << "Not defined" << std::endl;
	exit(-1);
      }
      else if ( param.Rotate ){
	if ( param.Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
	  sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
	else {
	  for ( j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
	objvar.mutateRotate( sigma );
      }
      else if ( param.GSA ){
	if ( param.Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
	  sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
	else {
	  for ( j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
	objvar.mutateGSA( sigma );
      }
      else if ( param.IDA ){
	if ( param.Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
	  sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
	else {
	  for ( j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( j = 0; j < sigma.size( ); j++ ){
	    sigma[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 1 ] )[ j ];
	  }
	}
	objvar.mutateIDA( sigma );
      }
      else if ( param.CMA ){
	if ( param.Recombine ){
	  objvar.recombineDiscrete       ( pap[ 0 ], mom[ 0 ] );
	  sigma .recombineGenIntermediate( pap[ 1 ], mom[ 1 ] );
	}
	else {
	  for ( j = 0; j < objvar.size( ); j++ ){
	    objvar[ j ] = dynamic_cast< ChromosomeT< double >& >( pap[ 0 ] )[ j ];
	  }
	  for ( j = 0; j < sigma.size( ); j++ ){
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
    for ( j = 0; j < param.LAMBDA; j++ ){
      switch(param.MOOF)
	{
	case 1:
	  f1 = SphereF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  f2 = SphereF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  break;
	case 2:
	  f1 = DebConvexF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  f2 = DebConvexF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  break;
	case 3:
	  f1 = DebConcaveF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  f2 = DebConcaveF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  break;
	case 4:
	  f1 = DebDiscreteF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  f2 = DebDiscreteF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  break;
	case 5:
	  f1 = FonsecaConcaveF1( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  f2 = FonsecaConcaveF2( dynamic_cast< ChromosomeT< double >& >( offsprings[j][0] ));
	  break;
	default:
	  std::cerr << "Specify fitness function. Current setting of MOOF (1-5):" << param.MOOF << std::endl;
	  break;
	}
      offsprings[j].setMOOFitnessValues( f1, f2 );
      XO[ j ] = f1;
      YO[ j ] = f2;
    }
    offsprings.aggregation( weight );
   
    // Selection ( Mu, Labmda ) or ( Mu + Lambda )
    parents.selectMuLambda( offsprings, NELITE );


    for ( i = 0; i < param.MU; i++ ){
      switch(param.MOOF)
	{
	case 1:
	  f1 = SphereF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = SphereF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 2:
	  f1 = DebConvexF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = DebConvexF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 3:
	  f1 = DebConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = DebConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 4:
	  f1 = DebDiscreteF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = DebDiscreteF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	case 5:
	  f1 = FonsecaConcaveF1( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  f2 = FonsecaConcaveF2( dynamic_cast< ChromosomeT< double >& >( parents[i][0] ));
	  break;
	default:
	  std::cerr << "Specify fitness function. Current setting of MOOF (1-5):" << param.MOOF << std::endl;
	  break;
	}
      parents[i].setMOOFitnessValues( f1, f2 );
      XP[ i ] = f1;
      YP[ i ] = f2;
    }
    parents.aggregation( weight );

    // Store Parent Data
    for ( j = 0; j < parents.size( ); j++ ){
      XP[ j ] = parents[ j ].getMOOFitness( 0 );
      YP[ j ] = parents[ j ].getMOOFitness( 1 );
    }

    // Archive
    for ( j = 0; j < offsprings.size( ); j++ ){
      domination = archive.Dominate( offsprings[ j ] );
      mindistanceIA = archive.distanceOnFitness( offsprings[ j ] );
      if ( domination > 2 ){
	archive.delDominateArchive( offsprings[ j ] );
	archive.addArchive( offsprings[ j ] );
      }
      else if ( domination == 2 ){
	if ( mindistanceIA > param.NICHE ){
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
    for ( j = 0; j < archive.size( ); j++ ){
      XA[ j ] = archive.readArchive( j ).getMOOFitness( 0 );
      YA[ j ] = archive.readArchive( j ).getMOOFitness( 1 );
    }

    if ( t % param.Interval == 0 ){
      if ( param.Stopping ){
	std::cout << "\n Hit Return Key for Continue" << std::endl;
	fgetc(stdin);
      }
    }
  }
  
  char dummy[ 80 ];
  unsigned filenumber = param.FileName;
  sprintf( dummy, "Archive%i.txt", filenumber );
  archive.saveArchive( dummy );
  std::cout << "\nThe result is stored in " << dummy << "\n" << std::endl;
  std::cout << "Number of Solutions" << std::endl;
  std::cout << "Number of Objectives" << std::endl;
  std::cout << "Solution 1 (f1,f2,...)" << std::endl;
  std::cout << "Solution 2 (f1,f2,...)\n" << std::endl;

  std::cout << "Hit Return Key for End" << std::endl;
  fgetc(stdin);

  delete [] XP;
  delete [] YP;
  delete [] XO;
  delete [] YO;
  delete [] XA;
  delete [] YA;

  return 0;
}


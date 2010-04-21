/*! cmaplus.cpp
 * ======================================================================
 *
 *  \file    :  cmaplus.cpp
 *  \date    :  15.12.2004
 *
 *  \author Copyright (c) 2005 Stefan Roth
 *
 *  Institut fuer Neuroinformatik
 *  Ruhr-Universitaet Bochum
 *  44780 Bochum, Germany<BR>
 *  Phone: +49-234-32-25558<BR>
 *  Fax:   +49-234-32-14209<BR>
 *  eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *  www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *
 *  \par Project:
 *      EALib
 *  \par Language and Compiler:
 *     C++, egcs (Linux)
 *
 *  \par File and Revision:
 *       $RCSfile: cmaplus.cpp,v $<BR>
 *
 *  \par Changes:
 *        $Log: cmaplus.cpp,v $
 *        Revision 2.8  2006/04/30 21:09:40  christian_igel
 *        removed show commands
 *
 *        Revision 2.7  2005/10/11 12:01:30  christian_igel
 *        simplified
 *
 *        Revision 2.6  2005/05/23 18:10:03  christian_igel
 *        changes for gcc 3.4
 *
 *        Revision 2.5  2004/12/30 14:17:03  shark-admin
 *          comment added to MyParam::readParam()
 *
 *        Revision 2.4  2004/12/16 15:52:04  shark-admin
 *        comments removed
 *
 *        Revision 2.3  2004/12/16 11:29:58  shark-admin
 *        formated
 *
 *        Revision 2.2  2004/12/16 11:28:11  shark-admin
 *        parameter class added
 *
 *        Revision 2.1  2004/12/15 17:13:09  shark-admin
 *        example 'cmaplus' entered repository: (1+lambda) CMA
 *
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
#include <EALib/ChromosomeCMA.h>
#include <ReClaM/Params.h>
#include <ReClaM/Paraboloid.h>

using namespace std;


// My own derived class for managing my configuration files:
//
class MyParams : public Params 
{
public:
  MyParams( int argc, char **argv ) : Params( argc, argv ) { 
    setDefault();
    if( confFile.length( ) > 0 ) scanFrom( confFile );
  }
	
  ~MyParams( ) { }
	
  void io( istream& is, ostream& os, FileUtil::iotype type )  {
    FileUtil::io( is, os, "n"     , n     , 5u       , type);
    FileUtil::io( is, os, "lambda", lambda, 1u       , type);
    FileUtil::io( is, os, "T"     , T     , 3500u    , type);
    FileUtil::io( is, os, "sigma" , sigma , 1.       , type);
    FileUtil::io( is, os, "a"     , a     , 1000.    , type);    
    FileUtil::io( is, os, "min"   , min   , -1.      , type);
    FileUtil::io( is, os, "max"   , max   , 1.       , type);
    FileUtil::io( is, os, "seed"  , seed  , 10u      , type);

  }
	
  // Method to show the current content of the class variable: 
  void monitor( ) {
    cout << "n      (problem dimension)            " << n	    << endl;
    cout << "lambda (offspring pop size)           " << lambda << endl;
    cout << "T      (max iterations)               " << T      << endl;
    cout << "sigma  (initial global step size)     " << sigma  << endl;
    cout << "a      (paraboloidal scale parameter) " << a	    << endl;
    cout << "min    (min initial allele value)     " << min    << endl;
    cout << "max    (max initial allele value)     " << max    << endl;
    cout << "seed   (of random number generator)   " << seed   << endl;
  }

  unsigned n;
  unsigned lambda;
  unsigned T;
  unsigned seed;
  double sigma;
  double a;
  double min;
  double max;
};


int main(int argc, char* argv[])
{
  MyParams param( argc, argv );
	
  int lambdaSucc;		  	
  unsigned i, t;
  Array<double> B;


  Rng::seed(param.seed);

  // init fitness function with random coordinate system
  Paraboloid f(param.n, param.a, true);

  ChromosomeCMA chrome(param.n);
  Individual    parentIndiv(chrome);
  Population    parent(1, parentIndiv), offspring(param.lambda, parentIndiv);
	
  parent.setMinimize();
  offspring.setMinimize();

  // init single parent
  (dynamic_cast<ChromosomeCMA&>(parent[0][0])).init(param.n, param.sigma, param.min,param.max, param.lambda);
  parent[0].setFitness(f.error(dynamic_cast<ChromosomeT<double>&>(parent[0][0])));

  for(t=0; t<param.T; t++) {
    lambdaSucc = 0; 
    parentIndiv = parent[0];
    for(i=0; i<param.lambda; i++) {
      offspring[i] = parentIndiv;
      (dynamic_cast<ChromosomeCMA&>(offspring[i][0])).mutate();
      offspring[i].setFitness(f.error(dynamic_cast<ChromosomeT<double>&>(offspring[i][0])));
      if(offspring[i].getFitness() < parentIndiv.getFitness()) {
	lambdaSucc++;
	if(offspring[i].getFitness() < parent[0].getFitness()) parent[0] = offspring[i];
      }
    }

    cout << t << " " 	<< parentIndiv.getFitness() << " " << lambdaSucc << " " << endl;
    
    if(lambdaSucc)  (dynamic_cast<ChromosomeCMA&>(parent[0][0])).updateCovariance(parentIndiv);
    (dynamic_cast<ChromosomeCMA&>(parent[0][0])).updateGlobalStepsize(lambdaSucc);
  }

#define VERBOSE
#ifdef VERBOSE
  cout << endl;
  param.monitor();
  ChromosomeCMA& solution = (dynamic_cast<ChromosomeCMA&>(parent[0][0]));
  cout << "\n\nFinal parametervector" << endl;
  for(i=0; i< parent[0][0].size(); i++) cout << solution[i] << " ";
#endif

  return EXIT_SUCCESS;
}


/*! cmaExample.cpp
* ======================================================================
*
*  \file    :  cmaExample.cpp
*  \date    :  26.05.2004
*
*  \author Copyright (c) 2004 Christian Igel
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
*       $RCSfile: cmaExample.cpp,v $<BR>
*
*  \par Changes:
*        $Log: cmaExample.cpp,v $
*        Revision 2.6  2005/08/03 12:00:14  christian_igel
*        changes limits/values includes etc.
*
*        Revision 2.5  2004/12/10 17:38:59  christian_igel
*        Removed unnecessary loop.
*
*        Revision 2.4  2004/06/17 23:05:26  saviapbe
*        The standard file header was for doxygen adapted.
*
*        Revision 2.3  2004/06/09 17:03:26  saviapbe
*        "cmaExmaple" was for VC++ 6.0 and VC++ v7.0 adapted.
*
*        Revision 2.2  2004/06/01 12:18:15  saviapbe
*
*        "#include <limits>" was changed to "#include <limits.h>".
*
*        Revision 2.1  2004/05/27 15:10:04  saviapbe
*
*        An example for CMA2004 was added.
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

namespace std {} using namespace std;
#include <EALib/cma2004.h>
#include <ReClaM/Paraboloid.h>
#include <iomanip>
#include <Array/ArrayOp.h>

#ifdef __SOLARIS__
#include <values.h>
#else
#include <limits>
#endif

#include <EALib/sqr.h>
double sphere(const vector<double> &v, double cond) {
  double sum = 0.;
  for(unsigned i = 0; i < v.size(); i++) 
    sum += sqr(pow(cond, double(i) / double(v.size() - 1)) * v[i]);
  return sum;
} 

bool isFeasible(const vector<double> &c) {
  return true;
}

int main() {
  unsigned i, j;
  //
  // fitness function
  //
  const unsigned Dimension  = 100;
  const unsigned a          = 100;  // determines problem condition
  const bool     rotate     = true;
  Paraboloid parabel(Dimension, a);
  parabel.init(0, 0., 0., rotate);

  //
  // EA parameters
  //
  const unsigned Iterations     = 1000;
  const double   MinInit        = -3.;
  const double   MaxInit        = 7.;
  const double   GlobalStepInit = 1.;
#ifndef __SOLARIS__
  const double   Worst          = std::numeric_limits< double >::max( );
#else
  const double   Worst          = MAXDOUBLE;
#endif

  //
  // strategy parameters
  //
  CMA2004 cma;

  unsigned Lambda             = cma.suggestLambda(Dimension);
  unsigned Mu                 = cma.suggestMu(Lambda);

  unsigned feasible;

  cerr << "(" << Mu << "/" << Mu << ", " << Lambda << ")-CMA-ES on " 
       << Dimension << "-dimensional paraboloid" << endl;
    
  //
  // initialize random number generator
  //
  Rng::seed(10);

  //
  // define parent populations
  //
  Population parents (Mu,     
                      ChromosomeT< double >( Dimension ),
                      ChromosomeT< double >( Dimension ));

  //
  // define offspring population
  //
  Population offsprings (Lambda,     
                         ChromosomeT< double >( Dimension ),
                         ChromosomeT< double >( Dimension ));

  //
  // we want to minimize the target function
  //
  offsprings.setMinimize( );
  parents   .setMinimize( );  

  //
  // initialize parent populations center of gravity
  //
  cerr << "init in CMA" << endl;

  static_cast< ChromosomeT< double >& >( parents[0][0] ).initialize(MinInit, MaxInit);

  for(j = 1; j < parents.size( ); ++j ) 
    static_cast< ChromosomeT< double >& >(parents[j][0]) =
      static_cast< ChromosomeT< double >& >(parents[0][0]);

  bool stop = false;
  cout << 0 << " "  
       << sphere(dynamic_cast< vector< double >& >(parents[0][0]), a) << endl;


  //
  // strategy parameters
  //

  vector< double > variance( Dimension );
  for(i = 0; i < Dimension; i++) variance[i] = 1.;

  cma.init(Dimension, variance, GlobalStepInit, parents, CMA2004::superlinear, CMA2004::rankmu); 
  //cma.init(Dimension, variance, GlobalStepInit, parents, CMA::linear, CMA::rankmu); 

  // the previous three lines are equal to 
  //
  // cma.init(Dimension, 1., parents);
  //
  // but demonstrate how the CMA can be initialized with different variances

  //
  // iterate
  //
  for( unsigned t = 0; t < Iterations; ++t ) {
    //
    // generate at least Mu + 1 feasible new offsprings
    //
    feasible = 0;
    for(unsigned k = 0; k < offsprings.size( ); k++ ) {
      do {
	cma.create(offsprings[k]);
	if(isFeasible(dynamic_cast< vector< double >& >(offsprings[k][0]))) {
	  feasible++;
	  offsprings[k].setFitness(sphere(dynamic_cast< vector< double >& >(offsprings[k][0]), a));
	  break;
	} else {
	  offsprings[k].setFitness(Worst);
	}
      } while(feasible < Mu + 1);
    }
    
    //
    // select (mu,lambda) or (mu+lambda)
    //
    parents.selectMuLambda(offsprings, 0u);
  
    //
    // update strategy parameters
    //
    cma.updateStrategyParameters(parents);

    if((!stop) && (parents.best( ).fitnessValue( ) < 10E-10)) {
      stop = true;
      cerr << "stop: " << (t + 1) * Lambda << endl;
    }
      
    cout << setprecision(14) <<  (t + 1) * Lambda << " "  << parents.best( ).fitnessValue( ) << endl;//" " << cma.getCVar() << endl;
  }
  return 0;
}

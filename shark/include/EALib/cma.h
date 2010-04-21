//===========================================================================
/*!
 *  \file cma.h
 *
 *
 *  \par Copyright (c) 1998-2003:
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
 *      $RCSfile: cma.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: cma.h,v $
 *      Revision 2.2  2005/10/06 09:25:49  christian_igel
 *      sigma -> delta
 *
 *      Revision 2.1  2005/01/25 17:23:50  shark-admin
 *      lower bound on delta * min(lambda): fabs(lambda(Dimension - 1)) -> sqrt(fabs(lambda(Dimension - 1)))
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




//=============================================================
//
// CMA stuff
//

#include <EALib/sqr.h>
//#include <MathConst.h>
#include <EALib/Population.h>
#include <Array/Array.h>
#include <Array/ArrayIo.h>
#include <LinAlg/linalg.h>

//
// initialize s, sDelta, C
//
class CMA {
public:
  ~CMA() {}
  CMA()  {}

  void estimateChi_n(unsigned a = 10000) {
    chi_n = 0;
    for(i = 0; i < a; i++) {
      double sum = 0;
      for(j = 0; j < Dimension; j++)
        sum += sqr(Rng::gauss(0,1));
      chi_n += sqrt(sum);
    }
    chi_n /=a;
  }

  void init(unsigned _Dimension, std::vector<double > var, double _delta) {
    Dimension = _Dimension;
    c        = 1. / sqrt( double( Dimension ) );
    cu       = sqrt((2.-c) * c);
    ccov     = 2. / (sqr(double( Dimension) ) + Dimension);
    D        = 1 / sqrt(double( Dimension ));
    chi_n    = sqrt( double(Dimension) ) * (1 - 1. / (4.*Dimension) +  1. /
(21.*sqr( double(Dimension))));

    z.resize         (Dimension);
    s.resize         (Dimension);
    sDelta.resize    (Dimension);
    C.resize         (Dimension, Dimension);
    B.resize         (Dimension, Dimension);
    lambda.resize    (Dimension);

    theMatrix.resize (Dimension, Dimension);
    theVector.resize (Dimension);

    delta = _delta;
    for (i = 0; i < Dimension; i++) {
      s(i) = sDelta(i) = 0.;
      for (j = 0; j < Dimension; j++) {
        if(i!=j) C(i,j) = 0;
        else C(i,j) = var[i];
      }
    }

    //
    // eigenvalues lambda and eigenvector matrix B
    //
    eigensymm(C, B, lambda);
  }

  void init(unsigned _Dimension, double _delta) {
    Dimension = _Dimension;
    c        = 1. / sqrt( double( Dimension ) );
    cu       = sqrt((2.-c) * c);
    ccov     = 2. / (sqr(double(Dimension)) + Dimension);
    D        = 1 / sqrt(double(Dimension));
    chi_n    = sqrt(double(Dimension)) * (1 - 1. / (4.*Dimension)
                                  +  1. / (21.*sqr(double(Dimension))));

    z.resize         (Dimension);
    s.resize         (Dimension);
    sDelta.resize    (Dimension);
    C.resize         (Dimension, Dimension);
    B.resize         (Dimension, Dimension);
    lambda.resize    (Dimension);

    theMatrix.resize (Dimension, Dimension);
    theVector.resize (Dimension);

    delta = _delta;
    for (i = 0; i < Dimension; i++) {
      s(i) = sDelta(i) = 0.;
      for (j = 0; j < Dimension; j++) {
        if(i!=j) C(i,j) = 0;
        else C(i,j) = 1;
      }
    }

    //
    // eigenvalues lambda and eigenvector matrix B
    //
    eigensymm(C, B, lambda);
  }

  //
  // calculate center
  //
  void center(ChromosomeT<double >& x, Population &p, unsigned c = 0) {
    unsigned i, j;
    unsigned Dimension = dynamic_cast<ChromosomeT<double >&
>(p[0][c]).size();
    for(j = 0; j < Dimension; j++) x[j] = dynamic_cast<ChromosomeT<double >&
>(p[0][c])[j] ;
    for(i = 1; i < p.size(); i++) {
      for(j = 0; j < Dimension; j++) {
        x[j] += dynamic_cast< ChromosomeT< double >& >(p[i][c])[j];
      }
    }
    for(j = 0; j < Dimension; j++) x[j] /= p.size();
  }

  void mutate(ChromosomeT<double> &x, Individual &o) {
    unsigned i, j;

    // draw random vector
    for(i = 0; i < Dimension; i++) {
      z(i) = Rng::gauss(0,1);
      dynamic_cast< ChromosomeT< double >& >(o[1])[i] = z(i);
    }

    // mutate objective variables, Eq. (1)
    for (i = 0; i < Dimension; i++) dynamic_cast< ChromosomeT< double >&
>(o[0])[i] = x[i];

    for (i = 0; i < Dimension; i++)
      for (j = 0; j < Dimension; j++) {
        dynamic_cast< ChromosomeT< double >& >(o[0])[i] += delta * B(i,j) *
sqrt(fabs(lambda(j))) * z(j);
      }
  }

  void mutate97(ChromosomeT<double> &x, Individual &o) {
    unsigned i, j;

    // draw random vector
    for(i = 0; i < Dimension; i++) {
      z(i) = Rng::gauss(0,1);
    }

    // mutate objective variables, Eq. (1)
    for (i = 0; i < Dimension; i++) dynamic_cast< ChromosomeT< double >&
>(o[0])[i] = x[i];

    for (i = 0; i < Dimension; i++)
      for (j = 0; j < Dimension; j++)
        dynamic_cast< ChromosomeT< double >& >(o[0])[i] += delta * B(i,j) *
sqrt(fabs(lambda(j))) * z(j);
  }

  void updateStrategyParameters(ChromosomeT<double> &x, Population &p,
                                double lowerBound = .0) {
    double normS;
    ChromosomeT<double> xPrime(Dimension);
    center(xPrime, p);

    // Eq. (2)
    for (i = 0; i < Dimension; i++) {
      s(i) = (1 - c) * s(i) + cu * sqrt(double(p.size())) / delta *
(xPrime[i] - x[i]);
    }
    // Eq. (3)
    for (i = 0; i < Dimension; i++)
      for (j = 0; j < Dimension; j++)
        C(i, j ) = ( 1 - ccov ) * C(i, j) + ccov * s(i) * s(j);


    // Eq. (4)
    normS     = 0;
    theVector = 0;

    Array<double> help(Dimension);
    help = 0;
    for (i = 0; i < Dimension; i++) {
      for (j = 0; j < p.size(); j++) {
        help(i) += dynamic_cast< ChromosomeT< double >& >(p[j][1])[i];
      }
    }

    for (i = 0; i < Dimension; i++)
      for (j = 0; j < Dimension; j++)
        theVector(i) += B(i,j) * help(j) / sqrt(double(p.size()));

    for (i = 0; i < Dimension; i++) {
      sDelta(i) = (1 - c) * sDelta(i) + cu * theVector(i);
      normS += sqr(sDelta(i));
    }
    normS  = sqrt( normS );

    // Eq. (5)
    delta *= exp(D * ( normS - chi_n ) / chi_n );

    //
    // eigenvalues lambda and eigenvector matrix B
    //
    eigensymm(C, B, lambda);

    //
    // lower bound
    //
    if ((delta * sqrt(fabs(lambda(Dimension - 1)))) < lowerBound) {
      delta = lowerBound / sqrt(fabs(lambda(Dimension - 1)));
    }
  }

  void updateStrategyParameters97(ChromosomeT<double> &x, Population &p) {
    double normS;
    ChromosomeT<double> xPrime(Dimension);
    center(xPrime, p);

    // Eq. (2)
    for (i = 0; i < Dimension; i++)
      s(i) = (1 - c) * s(i) + cu * sqrt(double(p.size())) / delta *
	(xPrime[i] - x[i]);

    // Eq. (3)
    for (i = 0; i < Dimension; i++)
      for (j = 0; j < Dimension; j++)
        C(i, j ) = ( 1 - ccov ) * C(i, j) + ccov * s(i) * s(j);

    // Eq. (4)
    normS      = 0;
    theMatrix = 0;
    theVector = 0;

    for (i = 0; i < Dimension; i++)
      for (j = 0; j < Dimension; j++)
        for (k = 0; k < Dimension; k++)
          theMatrix(i,j) += B(i,k) * B(j,k) / sqrt(fabs(lambda(k)));
    for (i = 0; i < Dimension; i++)
      for (j = 0; j < Dimension; j++)
        theVector(i) += theMatrix(i, j) * (xPrime[j] - x[j]);

    for (i = 0; i < Dimension; i++) {
      sDelta(i) = (1 - c) * sDelta(i) + cu * theVector(i) *
sqrt(double(p.size())) / delta;
      normS += sqr(sDelta(i));
    }
    normS  = sqrt( normS );

    // Eq. (5)
    delta *= exp(D * ( normS - chi_n ) / chi_n );

    //
    // eigenvalues lambda and eigenvector matrix B
    //
    eigensymm(C, B, lambda);
  }

  double getDelta()                 { return delta; }
  void   setDelta(double x)         { delta = x; }
  void   cutDelta(double x = 1E-10) { if(delta < x) delta = x; }
  double getChi_n()                 { return chi_n; }
  void   setChi_n(double x)         { chi_n = x; }
  double getCondition()             { return lambda(0) / lambda(Dimension -
1); }
  const  Array<double> &getC()      { return C; }
  const  Array<double> &getLambda() { return lambda; }

private:
  unsigned Dimension;
  unsigned i, j, k;
  double delta;
  double chi_n;
  double c;
  double cu;
  double ccov;
  double D;

  Array<double> z;
  Array<double> s;
  Array<double> sDelta;
  Array<double> C;
  Array<double> lambda;
  Array<double> B;
  Array<double> theMatrix;
  Array<double> theVector;
};

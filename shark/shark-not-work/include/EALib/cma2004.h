//===========================================================================
/*!
 *  \file CMA2004.h
 *
 *
 *  \par Copyright (c) 1998-2004:
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
 *      Shark
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: cma2004.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: cma2004.h,v $
 *      Revision 2.3  2005/01/25 17:23:50  shark-admin
 *      lower bound on sigma * min(lambda): fabs(lambda(Dimension - 1)) -> sqrt(fabs(lambda(Dimension - 1)))
 *
 *      Revision 2.2  2004/12/22 16:02:43  christian_igel
 *      Corrected Bug in suggestMu.
 *
 *      Revision 2.1  2004/06/18 03:47:25  saviapbe
 *      The dummy linking to examples.
 *
 *      Revision 2.0  2004/05/26 19:53:16  saviapbe
 *
 *      Revision tag reset to revision tag 2.x
 *      The class CMA2004 was added to Shark.
 *
 *
 *  This file is part of Shark. This library is free software;
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

#include <EALib/sqr.h>
#include <EALib/Population.h>
#include <Array/ArrayOp.h>
#include <Array/ArrayIo.h>
#include <LinAlg/linalg.h>
/*!
*
* \example cmaExample.cpp
* 
*/
class CMA2004 {
public:
  ~CMA2004() {}
  CMA2004()  {}

  typedef enum { equal, linear, superlinear } RecombType;
  typedef enum { source, source2, paper }     InitType;
  typedef enum { rankone, rankmu }            UpdateType;

  unsigned suggestLambda(unsigned dimension) {
    unsigned lambda             = unsigned( 4. + floor(3. * log((double) dimension)) );
    // heuristics for small search spaces
    if(lambda > dimension) lambda = dimension; // CI's golden rule :-)
    if(lambda < 5) lambda = 5; // Hansen & Ostermeier's lower bound
    return lambda;
  }

  unsigned suggestMu(unsigned lambda, RecombType recomb = superlinear) { 
    if(recomb == equal) return  unsigned( floor(lambda / 4.) ); 
    return  unsigned( floor(lambda / 2.) );
  }
    
  void init(unsigned dimension, 
	    std::vector<double > var, double _sigma, 
	    Population &p, 
	    RecombType recomb  = superlinear, 
	    UpdateType cupdate = rankmu, 
	    InitType init      = paper) {

    unsigned mu = p.size();

    n     = dimension;
    sigma = _sigma;

    w.resize         (mu);
    x.resize         (n);
    xPrime.resize    (n);
    z.resize         (n);
    pc.resize        (n);
    ps.resize        (n);
    Z.resize         (n, n);
    C.resize         (n, n);
    B.resize         (n, n);
    lambda.resize    (n);
    theVector.resize (n);
    meanz.resize     (n);

    switch(recomb) {
    case equal:
      for(i = 0; i < mu; i++) w(i) = 1;
      break;
    case linear: 
      for(i = 0; i < mu; i++) w(i) = mu - i;
      break;
    case superlinear: 
      for(i = 0; i < mu; i++) w(i) = log(mu + 1.) - log(1. + i);
      break;
    }

    double wSum    = 0;
    double wSumSqr = 0;
    for(i = 0; i < mu; i++) {
      wSum += w(i);
      wSumSqr += sqr(w(i));
    }
    w /= wSum; // normalizing weights
    wSumSqr /= sqr(wSum);

    cc      = 4. / (4. + n);
    cs      = 10. / (20. + n);
    ccu     = sqrt((2.-cc) * cc);
    csu     = sqrt((2.-cs) * cs);
    chi_n   = sqrt( double(n) ) * (1 - 1. / (4. * n) +  1. / (21. * sqr( double(n))));
      
    mueff   = 1 / wSumSqr;
    mucov   = mueff;
    if(cupdate == rankone) mucov = 1.;
    d       = max(1., 3*mueff / (n + 10.)) + cs; // 1./cs;

    switch(init) {
    case source:
      ccov     = 1./mucov * 2./sqr(n + (2.)) 
	+ (1 - 1./mucov) * (2*mueff - 1) / (sqr(n + 2) + 2 * mueff);
      break;
    case paper:
      ccov     = 1./mucov * 2./sqr(n + sqrt(2.)) 
	+ (1 - 1./mucov) * min(1., (2*mueff - 1) / (sqr(n + 2) + mueff));
      break;
    default:
      cerr << "unknow CMA initialization command" << endl;
      exit(EXIT_FAILURE);
    }
    								      
    // init COG								      
    cog(x, p);

    // init paths
    for(i = 0; i < n; i++) {
      pc(i) = ps(i) = 0.;
      for(j = 0; j < n; j++) {
        if(i!=j) C(i,j) = 0;
        else C(i,j) = var[i];
      }
    }
    
    // eigenvalues lambda and eigenvector matrix B
    eigensymm(C, B, lambda);
  }

  void init(unsigned dimension, double _sigma,  Population &p) {
    std::vector<double> var(dimension);
    for(i = 0; i < dimension; i++) var[i] = 1;
    init(dimension, var, _sigma, p);
  }

  //
  // calculate weighted mean
  //
  void cog(ChromosomeT<double >& a, Population &p, unsigned c = 0) {
    SIZE_CHECK( n == dynamic_cast<ChromosomeT<double >& >(p[0][c]).size() );
    for(j = 0; j < n; j++) {
      a[j] = dynamic_cast< ChromosomeT< double >& >(p[0][c])[j] * w(0);
      for(i = 1; i < p.size(); i++) {
        a[j] += dynamic_cast< ChromosomeT< double >& >(p[i][c])[j] * w(i);
      }
    }
  }

  //
  // mutation after global intermediate recombination
  //
  void create(Individual &o) {
    for(i = 0; i < n; i++) {
      // draw random vector
      z(i) = Rng::gauss(0,1);
      dynamic_cast< ChromosomeT< double >& >(o[1])[i] = z(i);
      // global intermediate recombination, Eq. (1)
      dynamic_cast< ChromosomeT< double >& >(o[0])[i] = x[i];
    }

    // mutate objective variables, Eq. (1)
    for(i = 0; i < n; i++)
      for(j = 0; j < n; j++) 
        dynamic_cast< ChromosomeT< double >& >(o[0])[i] += sigma * B(i,j) * sqrt(fabs(lambda(j))) * z(j);
  }

  //
  // do the CMA
  //
  void updateStrategyParameters(Population &p, double lowerBound = .0) {
    double normPS = 0.;

    // COG of new parents
    cog(xPrime, p);
    cog(meanz , p, 1);

    theVector = 0;
    for(i = 0; i < n; i++)
      for(j = 0; j < n; j++)
        theVector(i) += B(i,j) * meanz[j];

    // Eq. (2) & Eq. (4)
    for(i = 0; i < n; i++) {
      pc(i) = (1 - cc) * pc(i) + ccu * sqrt(mueff) / sigma * (xPrime[i] - x[i]);
      ps(i) = (1 - cs) * ps(i) + csu * sqrt(mueff) * theVector(i);
      normPS += sqr(ps(i));
    }
    normPS  = sqrt( normPS );

    // Eq. (3)
    for(i = 0; i < n; i++) 
      for(j = 0; j < n; j++)
	for(Z[i][j] = 0., k = 0; k < p.size(); k++)
	  Z(i, j ) += w(k) * (dynamic_cast< ChromosomeT< double >& >(p[k][0])[i] - x[i]) *
	    (dynamic_cast< ChromosomeT< double >& >(p[k][0])[j] - x[j]);

    for(i = 0; i < n; i++) 
      for(j = 0; j < n; j++)
        C(i, j ) = (1 - ccov) * C(i, j) + ccov * 
	  (1./mucov * pc(i) * pc(j) + (1. - 1./mucov) * 1./sqr(sigma) * Z(i, j));

    // Eq. (5)
    sigma *= exp((cs / d) * (normPS / chi_n - 1));

    // eigenvalues lambda and eigenvector matrix B
    eigensymm(C, B, lambda);

    // lower bound
    if ((sigma * sqrt(fabs(lambda(n - 1)))) < lowerBound) 
      sigma = lowerBound / sqrt(fabs(lambda(n - 1)));

    // new COG becomes old COG
    x = xPrime;
  }

  double getSigma()                 { return sigma; }
  void   setSigma(double x)         { sigma = x; }
  void   setChi_n(double x)         { chi_n = x; }
  double getCondition()             { return lambda(0) / lambda(n-1); }
  const  Array<double> &getC()      { return C; }
  const  Array<double> &getLambda() { return lambda; }

private:
  unsigned n;
  unsigned i, j, k;
  double sigma;
  double chi_n;
  double cc;
  double cs;
  double csu;
  double ccu;
  double ccov;
  double d;
  double mueff;
  double mucov;

  ChromosomeT<double> x;
  ChromosomeT<double> xPrime;
  ChromosomeT<double> meanz;

  Array<double> z;
  Array<double> pc;
  Array<double> ps;
  Array<double> C;
  Array<double> Z;
  Array<double> lambda;
  Array<double> B;
  Array<double> w;
  Array<double> theVector;
};

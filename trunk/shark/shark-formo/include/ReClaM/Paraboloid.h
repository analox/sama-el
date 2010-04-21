//===========================================================================
/*!
 *  \file Paraboloid.h
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
 *      ReClaM
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: Paraboloid.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: Paraboloid.h,v $
 *      Revision 2.5  2006/04/20 13:28:59  glasmtbl
 *      conflicting versions 2.3 and 2.4 merged
 *
 *      Revision 2.3  2005/10/11 12:02:10  christian_igel
 *      simplified
 *
 *      Revision 2.1  2004/06/01 15:28:35  saviapbe
 *
 *      The names od parameters were changed.
 *
 *      Revision 2.0  2004/05/26 19:56:44  saviapbe
 *
 *      Revision tag reset to revision tag 2.x
 *      The class "Paraboloid" was added to Shark and is as from now available as a test function.
 *
 *
 *  This file is part of ReClaM. This library is free software;
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
// paraboloid test function
//

#ifndef PARABOLOID_H
#define PARABOLOID_H
#include <ReClaM/ModelInterface.h>
#include <Rng/GlobalRng.h>
#ifndef _WIN32
#include <values.h>
#endif
#include <math.h>
#include <EALib/sqr.h>
#include <vector>
#include <Array/ArrayOp.h>
#include <Array/ArrayIo.h>

class Paraboloid : virtual public ModelInterface {
 public:
  Paraboloid(unsigned _d, double _c = 10, bool basis = true) {
    unsigned i;
    d = _d;
    cond = _c;

    w   .resize(d);
    dmdw.resize(1, d);
    dedw.resize(d);
    B   .resize(d, d);
    w  = 0;
    dmdw = 0;
    dedw = 0;
    if(basis) generateBasis();
    else {
      B  = 0;
      for(i = 0; i < d; i++) B(i, i) = 1;
    }
  }

  void model(const Array<double> &input, Array<double> &output) {
  }

  void dmodel(const Array<double> &input) {
  }

  double error(const std::vector<double> &v) {
    unsigned i;
    if(v.size() != d) {
      std::cerr << "dimension mismatch" << std::endl;
      throw "dimension mismatch";
    }
    double sum = 0.;
    for(i = 0; i < d; i++) w(i) = v[i];
    for(i = 0; i < d; i++) {
      sum += sqr( pow(cond, (double(i) / double(d - 1))) * scalarProduct(w, B.col(i)) ) ;
      // sum += sqr(pow(cond, double(i) / double(d - 1)) * scalarProduct(w, B[i]));
    }
    return sum;
  }    

  double error(const Array<double> &input, const Array<double> &target) {
    double sum = 0.;
    unsigned i;
    for(i = 0; i < d; i++) 
      sum += sqr(pow(cond, double(i) / double(d - 1)) * scalarProduct(w, B.col(i)));
    return sum;
  }

  double derror(const Array<double> &input, const Array<double> &target, bool returnError = true) {
    double sum = 0.;
    unsigned k, i; 

    dedw = 0;
    for(k = 0; k < d; k++) {
      sum += sqr(pow(1, (double) k / (double) (d - 1)) * scalarProduct(w, B.col(k)));
      for(i = 0; i < d; i++) 
        dedw(k) += 2 * pow(cond, 2 * double(i) / double(d - 1)) * scalarProduct(w, B.col(i)) * B(k, i);
    }
    return sum;
  }

  void init(unsigned s, double a, double b, bool basis = true) {
    unsigned i;
    Rng::seed(s);
    if(basis) generateBasis();
    for(i = 0; i < d; i++) w(i) = Rng::uni(a, b);
  }

  Array<double> &getBasis() { return B; }

protected:
  void generateBasis() {
    unsigned i, j, c;
    Array<double> H;
    B.resize(d, d);
    H.resize(d, d);
    for(i = 0; i < d; i++) {
      for(c = 0; c < d; c++) {
        H(i, c) = Rng::gauss(0, 1);
        //H(i, c) = (c==i) ? 1 : 0;
      }
    }
    B = H;
    for(i = 0; i < d; i++) {
      for(j = 0; j < i; j++) 
        for(c = 0; c < d; c++) 
          B(i, c) -= scalarProduct(H[i], H[j]) * H(j, c) / scalarProduct(H[j], H[j]);
      H = B;
    }
    for(i = 0; i < d; i++) {
      double normB = sqrt(scalarProduct(B[i], B[i]));
      for(j = 0; j < d; j++)
	B(i, j) = B(i, j) / normB;
    }
  }
  
  Array<double > B;
  unsigned       d;
  double         cond;
};

#endif




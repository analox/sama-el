//===========================================================================
/*!
 *  \file SquaredError.h
 *
 *  \brief Calculates the sum-of-squares error. 
 *
 *  These methods can be used as error measure model for the learning
 *  process.  
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Copyright (c) 1999-2001:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *
 *  \par Project:
 *      ReClaM
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: SquaredError.h,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:22 $
 *
 *  \par Changes:
 *      $Log: SquaredError.h,v $
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:27:41  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2001/11/30 15:02:49  rudi
 *      de -> dedw
 *      dw -> dmdw
 *      fdmodel (.,.) -> dmodel
 *      return value for derror
 *      fderror removed
 *      doxygen comments added
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

#ifndef SQUARED_ERROR_H
#define SQUARED_ERROR_H

#include <cmath>
#include "ReClaM/ModelInterface.h"

//===========================================================================
/*!
 *  \brief Calculates the sum-of-squares error. 
 *
 *  These methods can be used as error measure model for the learning
 *  process.  
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class SquaredError : virtual public ModelInterface {
 public:

//===========================================================================
/*!
 *  \brief Calculates the squared error between the output and the target
 *         vector.
 *
 *  Measures the euklidian distance between the model output \em model(in),
 *  calculated from the input vector \em in, and the target vector \em out. 
 *  Consider the case of a N-dimensional output vector, i.e. a neural network 
 *  with \em N output neurons, and a set of \em P patterns. In this case the 
 *  function calculates 
 *  \f[E = \sum_{p=1}^P\sum_{i=1}^N(model(in)_{ip} - out_{ip})^{2}\f]
 *
 *      \param  in Input vector for the model.
 *      \param  out Target vector.
 *      \return The squared error \em E.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
  double error(const Array<double> &in, const Array<double> &out) {
    double se = 0;
    if(in.ndim() == 1) {
      Array<double> output(out.dim(0));
      model(in, output);
      for(unsigned c = 0; c < out.dim(0); c++) {
	se += (out(c) - output(c)) * (out(c) - output(c));
      }
    } else {
      Array<double> output(out.dim(1));
      for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	model(in[pattern], output);
	for(unsigned c = 0; c < out.dim(1); c++) {
	  se += (out(pattern, c) - output(c)) * (out(pattern, c) - output(c));
	}
      }
    }
    return se;
  }

//===========================================================================
/*!
 *  \brief Calculates the derivatives of the squared error with respect
 *         to the parameters ModelInterface::w.
 *
 *  According to the equation in the description for the function
 *  #error the derivatives of the squared error can be calculated
 *  with
 *  \f[
 *  \frac{E}{w_j} = 2 \sum_{p=1}^P \sum_{i=1}^N (model(in)_{ip} - out_{ip})
 *  \frac{model(in)_{ip}}{w_j}
 *  \f]
 *  The results are written to the vector ModelInterface::dedw.
 *
 *  Usually, as a byproduct of the calculation of the derivative one gets the 
 *  the error \f$E\f$ itself very efficiently. Therefore, the method #error
 *  gives back this value. This additional effect can be switched of by means
 *  of the third parameter (returnError = false).
 *
 *      \param  in Input vector for the model.
 *      \param  out The target vector.
 *      \param  returnError Determines whether or not to calculate the error 
 *                          itself. By default the error is calculated.  
 *      \return The error \em E if \em returnError is set to "true", "-1"
 *              otherwise.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */

  double derror(const Array<double> &in, const Array<double> &out, bool returnError = true) {
    dedw = 0;
    double se = 0;
    if(in.ndim() == 1) {
      Array<double> output(out.dim(0));
      dmodel(in, output);
      for(unsigned c = 0; c < output.nelem(); c++) {
	if (returnError)
	  se += (out(c) - output(c)) * (out(c) - output(c));
	for(unsigned i = 0; i < dedw.nelem(); i++)
	  dedw(i) -= (out(c) - output(c)) * dmdw(c, i);
      }
      for(unsigned i = 0; i < dedw.nelem(); i++) dedw(i) *= 2;
    } else {
      Array<double> output(out.dim(1));
      for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	dmodel(in[pattern], output);
	for(unsigned c = 0; c < output.nelem(); c++) {
	  if(returnError)
	    se += (out(pattern, c) - output(c)) * (out(pattern, c) - output(c));
	  for(unsigned i = 0; i < dedw.nelem(); i++)
	    dedw(i) -= (out(pattern, c) - output(c)) * dmdw(c, i);
	}
      }
      for(unsigned i = 0; i < dedw.nelem(); i++) dedw(i) *= 2;
    }
    if (!returnError)
      se = -1;
    return se;
  }
};

#endif







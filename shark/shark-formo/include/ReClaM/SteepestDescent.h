//===========================================================================
/*!
 *  \file SteepestDescent.h
 *
 *  \brief The simplest learning strategy.
 *
 *  The steepest descent strategy uses the calculated error derivatives
 *  to adjust the weights of the network. 
 *         
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
 *      $RCSfile: SteepestDescent.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2006/04/10 15:48:14 $
 *
 *  \par Changes:
 *      $Log: SteepestDescent.h,v $
 *      Revision 2.1  2006/04/10 15:48:14  glasmtbl
 *      Error in the documentation corrected.
 *
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:27:41  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2001/11/30 15:03:51  rudi
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

#ifndef STEEPEST_DESCENT_H
#define STEEPEST_DESCENT_H

#include "Array/Array.h"
#include "Array/ArrayOp.h"
#include "ReClaM/ModelInterface.h"

//===========================================================================
/*!
 *  \brief The simplest learning strategy.
 * 
 *  The steepest descent strategy uses the calculated error derivatives
 *  to adjust the weights of the network. 
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
class SteepestDescent : virtual public ModelInterface {
 public:

//===========================================================================
/*!
 *  \brief Iterative updating of the weights by presenting all patterns.
 *
 *  Similar to #learnBatch, but here all patterns are presented to the
 *  network sequentially and the weights are updated after each pattern.
 *
 *      \param  in Input vector for the model.
 *      \param  target Target vector. 
 *      \param  lr Learning rate \f$ \eta \f$.
 *      \param  mu Momentum parameter \f$ \mu \f$.
 *      \return None.
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
  void learnBatch(const Array<double> &in, const Array<double> &target, double lr = 0.1, double mu = 0.3) {
    if(!m.ndim()) {
      m.resize(dedw);
      m = 0;
    }
    derror(in, target,false);
    m = -lr * dedw + mu * m;
    w += m;
  }

//===========================================================================
/*!
 *  \brief Iterative updating of the weights by presenting one pattern.
 *
 *  After presenting one pattern to the network the weights are iteratively
 *  updated by moving a short distance in the direction of the
 *  greatest decrease of the error, i.e. in the direction of the
 *  negative gradient by using the error derivatives #dedw.
 *  Additionally a momentum term is added to the gradient descent
 *  to deal with the problem of widely differing eigenvalues.
 *  The update of a weight at time \f$t\f$ to a weight at time
 *  \f$t+1\f$ is then given as
 *  \f[
 *      \omega_{t+1} = \omega_t + \eta \nabla f(\omega_t) + \mu \Delta w_t
 *  \f]
 *
 *      \param  in Input vector for the model.
 *      \param  target Target vector. 
 *      \param  lr Learning rate \f$ \eta \f$.
 *      \param  mu Momentum parameter \f$ \mu \f$.
 *      \return None.
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
  void learnOnline(const Array<double> &in, const Array<double> &target, double lr = 0.1, double mu = 0.3) {
    if(!m.ndim()) {
      m.resize(dedw);
      m = 0;
    }
    if(in.ndim() == 1) {
      derror(in, target);
      m = -lr * dedw + mu * m;
      w += m;
    } else {
      for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) { 
	derror(in[pattern], target[pattern],false);
	m = -lr * dedw + mu * m;
	w += m;
      }
    }
  }

  //===========================================================================
/*!
 *  \brief Clears the array of momentum values.
 *
 *  The momentum values are the values that are added to the weights
 *  for updating.
 *
 *  \return None.
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
  void clearMemory() { m = 0; };

 private:

  Array<double> m;

};

#endif



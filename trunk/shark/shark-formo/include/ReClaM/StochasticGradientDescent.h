//===========================================================================
/*!
 *  \file StochasticGradientDescent.h
 *
 *  \brief Learning strategies based on SteepestDescent, but with
 *         randomly chosen patterns.
 *
 *  Similar to the most simple learning strategy SteepestDescent,
 *  the weights of the network are iteratively updated to minimize
 *  the error. When working with more than one pattern and gradient
 *  descent, all patterns are sequentially presented to the network 
 *  and the weight update takes place after every pattern.
 *  Here you will find two methods that work on the base of
 *  SteepestDescent, but that will chose the next pattern for
 *  presentation by random. 
 *         
 *
 *  \author  P. Stagge
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
 *      $RCSfile: StochasticGradientDescent.h,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:22 $
 *
 *  \par Changes:
 *      $Log: StochasticGradientDescent.h,v $
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:27:41  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2001/11/30 15:05:26  rudi
 *      include-statements adapted to new version
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


#ifndef STOCHASTIC_GRADIENT_DESCENT_H
#define STOCHASTIC_GRADIENT_DESCENT_H

#include "Array/ArrayOp.h"
#include "ReClaM/ModelInterface.h"
#include "Rng/GlobalRng.h"

//===========================================================================
/*!
 *  \brief Learning strategies based on SteepestDescent, but with
 *         randomly chosen patterns.
 *
 *  Similar to the most simple learning strategy SteepestDescent,
 *  the weights of the network are iteratively updated to minimize
 *  the error. When working with more than one pattern and gradient
 *  descent, all patterns are sequentially presented to the network 
 *  and the weight update takes place after every pattern.
 *  Here you will find two methods that work on the base of
 *  SteepestDescent, but that will chose the next pattern for
 *  presentation by random. 
 *
 *  \author  P. Stagge
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class StochasticGradientDescent : virtual public ModelInterface {
public:

//===========================================================================
/*!
 *  \brief Iterative update of weights, patterns are randomly chosen,
 *         each pattern only once.
 *
 *  Similar to SteepestDescent, where the weights are iteratively
 *  updated by moving a short distance in the direction of the
 *  greatest decrease of the error, i.e. in the direction of the
 *  negative gradient by using the error derivatives #dedw.
 *  Additionally a momentum term is added to the gradient descent
 *  to deal with the problem of widely differing eigenvalues.
 *  The update of a weight at time \f$t\f$ to a weight at time
 *  \f$t+1\f$ is then given as
 *
 *  \f$
 *      \omega_{t+1} = \omega_t + \eta \nabla f(\omega_t) + \mu \Delta w_t
 *  \f$
 *
 *  In contrast to SteepestDescent the weights for more than
 *  one pattern are not updated by choosing the patterns sequentially,
 *  but randomly. Compared to #learnfast a pattern can not be chosen
 *  more than once.
 *
 *
 *      \param  in Input vector for the model.
 *      \param  target Target vector.
 *      \param  lr Learning rate \f$ \eta \f$.
 *      \param  mu Momentum parameter \f$ \mu \f$.
 *      \return None.
 *
 *  \author  P. Stagge
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
  void learn(const Array<double> &in, const Array<double> &target, double lr = 0.1, double mu = 0.3) {
    if(!m.ndim() || m.dim(0)!=dedw.dim(0)) {
    //if(!m.ndim()) {
    //geaendert, mt feb-2000
      m.resize(dedw);
      m = 0;
    }
    if(in.ndim() == 1) {
      derror(in, target,false);
      m = -lr * dedw + mu * m;
      w += m;
    } else {
      unsigned index, pattern;
      vector< int > vec( in.dim(0) );
      vector< int >::iterator it;
       for( unsigned i = 0; i < vec.size( ); ++i )
	vec[ i ] = i;
      for(unsigned number = 0; number < in.dim(0); ++number) { 
	index   = Rng::discrete(0,vec.size( )-1);
	pattern = vec[ index ];
	for(it = vec.begin(); *it != (int)pattern; it++);
	vec.erase(it);
	derror(in[pattern], target[pattern],false);
	m = -lr * dedw + mu * m;
	w += m;
      }
    }
  }


//===========================================================================
/*!
 *  \brief Iterative update of weights, patterns are randomly chosen.
 *
 *  Similar to #learn, but to save time it is not checked here whether a 
 *  pattern was already chosen.
 *
 *      \param  in Input vector for the model.
 *      \param  target Target vector.
 *      \param  lr Learning rate \f$ \eta \f$.
 *      \param  mu Momentum parameter \f$ \mu \f$.
 *      \return None.
 *
 *  \author  P. Stagge
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
  void learnfast(const Array<double> &in, const Array<double> &target, 
                 double lr = 0.1, double mu = 0.3                  ) 
  {
    if(!m.ndim()) {
      m.resize(dedw);
      m = 0;
    }
    if(in.ndim() == 1) {
      derror(in, target,false);
      m = -lr * dedw + mu * m;
      w += m;
    } else {
      unsigned pattern;
      for(unsigned number = 0; number < in.dim(0); ++number) { 
	pattern = Rng::discrete(0,in.dim(0)-1);
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
 *  \author  P. Stagge
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
  void clearMemory() 
  { 
      m = 0; 
  };


//===========================================================================
/*!
 *  \brief Clears the array of momentum values.
 *
 *  The momentum values are the values that are added to the weights
 *  for updating.
 *
 *  \return None.
 *
 *  \author  P. Stagge
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
  void clearMomentumMemory() 
  { 
      m = 0;
  };


 private:

  Array<double> m;  // Momentumterm

};

#endif



//===========================================================================
/*!
 *  \file Quickprop.h
 *
 *  \brief This file includes the popular heuristic optimization
 *         method named "Quickprop".
 *         
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Copyright (c) 1999-2000:
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
 *      $RCSfile: Quickprop.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2005/08/03 12:04:30 $
 *
 *  \par Changes:
 *      $Log: Quickprop.h,v $
 *      Revision 2.1  2005/08/03 12:04:30  christian_igel
 *      changes limits/values, MAXDOUBLE etc.
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
 *      Revision 1.3  2002/02/06 14:50:35  rudi
 *      Doxygen comments added.
 *
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

#ifndef QUICKPROP_H
#define QUICKPROP_H

#include "Array/Array.h"
#include "ReClaM/ModelInterface.h"
#ifndef __SOLARIS__
#include <limits>
#else
#include <values.h>
#endif

#include <fstream>
#include "EALib/sqr.h"
#include <math.h>


//===========================================================================
/*!
 *  \brief This class offers methods for using the popular heuristic
 *         "Quickprop" optimization algorithm. 
 *
 *  The "Quickprop" algorithm sees the weights of a network as if they were
 *  quasi-independant and tries to approximate the error surface, as a 
 *  function of each of the weights, by a quadratic polynomial (a parabola).
 *  Then two successive evaluations of the error function and an evaluation
 *  of its gradient follow to determine the coefficients of the
 *  polynomial. At the next step of the iteration, the weight parameter
 *  is moved to the minimum of the parabola. This leads to a calculation
 *  of the update factor for weight no. \f$i,\ i = 1, \dots, n\f$ from 
 *  iteration \f$t\f$ to iteration \f$t + 1\f$ as
 *
 *  \f$
 *  \Delta w_i^{(t + 1)} = \frac{g_i^{(t)}}{g_i^{t - 1} - g_i^{t}} 
 *                         \Delta w_i^{(t)}
 *  \f$
 *
 *  where \f$g_i^{(t)} = \frac{\partial E}{\partial w_i^{(t)}}\f$. <br>
 *
 *  For further information about this algorithm, please refer to: <br>
 *
 *  S. E. Fahlman, <br>
 *  "Faster-learning variations on back-propagation: an empirical study." <br>
 *  In D. Touretzky, G.E. Hinton and T,J, Sejnowski (Eds.), <br>
 *  "Proceedings of the 1988 Connectionist Models Summer School",
 *  pp. 38-51, San Mateo, CA: Morgan Kaufmann <br>
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 *
 */
class Quickprop : virtual public ModelInterface {
 public:


//===========================================================================
/*!
 *  \brief Prepares the Quickprop algorithm for the currently used network.
 *
 *  Internal variables of the class instance are initialized and memory
 *  for used structures adapted to the used network.
 *
 *  \return none
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
  void initQuickprop() 
  {
    deltaw.resize(w.nelem());
    dedwOld.resize(w.nelem());
    dedwOld = 0; deltaw = 0; 
  }


//===========================================================================
/*!
 *  \brief Performs one run of a modified Quickprop algorithm, that
 *         improves the performance.
 *
 *  The error \f$E^{(t)}\f$ of the used network for the current iteration 
 *  \f$t\f$ is calculated and the values of the weights \f$w_i,\ i = 1, 
 *  \dots, n\f$ are optimized depending on this error.
 *
 *  \param in    The input patterns used for the training.
 *  \param out   The target values for the input patterns.
 *  \param nu    The learning rate \f$\nu\f$, by default this is "1.5".
 *  \param dMax  Upper limit for the increase factor: No weight adaptation
 *               may be greater than the dMax times the previous 
 *               adaptation. The default value is "1.75".
 *  \return      Always "0" (not very reasonable).
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
  int quickprop(const Array<double> &in,
		const Array<double> &out,
		double              nu   = 1.5  ,   // learning rate
	        double              dMax = 1.75)     // upper limit of n
  {
    derror(in, out,false);

    for(unsigned  i = 0; i < w.nelem(); i++) {
      // avoid division by zero
      if(dedwOld(i) == dedw(i)) break;
      
      // quadratic part
      delta = dedw(i) / (dedwOld(i) - dedw(i)) * deltaw(i);;
      
      // add gradient term "unless the current slope is opposite in sign from
      // the previous slope"
      // initial step is steepest descent
      if(dedw(i) * dedwOld(i) >= 0) {
	delta -= nu * dedw(i);
      }
      
      // limit growth factor if slopes have the same sign
      if( (dedw(i) * dedwOld(i) > 0) &&
	  (fabs(delta) > fabs(dMax * deltaw(i))) )
	delta = sgn(delta) * fabs(dMax * deltaw(i));
      
      deltaw(i) = delta;
      
      w(i) += deltaw(i);
    }
    
    dedwOld = dedw;
    return 0;
  }


//===========================================================================
/*!
 *  \brief Performs one run of the original Quickprop algorithm.
 *
 *  The error \f$E^{(t)}\f$ of the used network for the current iteration 
 *  \f$t\f$ is calculated and the values of the weights \f$w_i,\ i = 1, 
 *  \dots, n\f$ are optimized depending on this error.
 *
 *  \param in    The input patterns used for the training.
 *  \param out   The target values for the input patterns.
 *  \param nu    The learning rate \f$\nu\f$, by default this is "1.5".
 *  \param dMax  Upper limit for the increase factor: No weight adaptation
 *               may be greater than the dMax times the previous 
 *               adaptation. The default value is "1.75".
 *  \return      Always "0" (not very reasonable).
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
  int quickpropOrig(const Array<double> &in,
		    const Array<double> &out,
		    double              nu   = 1.5  ,   // learning rate
		    double              dMax = 1.75)    // upper limit of n
  {
    derror(in, out,false);

    double shrinkFactor = dMax / (1 + dMax);

    for(unsigned  i = 0; i < w.nelem(); i++) {
      // avoid division by zero
      if(dedwOld(i) == dedw(i)) break;
      
      delta = 0.;
      
      if(deltaw(i) < 0.) {
	if(dedw(i) > 0.) 
	  // add negative gradient
	  delta -= nu * dedw(i);
	if(dedw(i) > shrinkFactor * dedwOld(i))
	  delta += dMax * deltaw(i);
	else
	  // quadratic part
	  delta += dedw(i) / (dedwOld(i) - dedw(i)) * deltaw(i);
      }
      else if(deltaw(i) > 0.) {
	if(dedw(i) < 0.) 
	  // add negative gradient
	  delta -= nu * dedw(i);
	if(dedw(i) < shrinkFactor * dedwOld(i))
	  delta += dMax * deltaw(i);
	else
	  // quadratic part
	  delta += dedw(i) / (dedwOld(i) - dedw(i)) * deltaw(i);
      } else {
	delta -= nu * dedw(i);
      }

      deltaw(i) = delta;
      w(i) += deltaw(i);
    }
    
    dedwOld = dedw;
    return 0;
  }


//===========================================================================
/*!
 *  \brief Resets internal variables used by the current class instance.
 *
 *  \return none
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
  void clearMemory() { dedwOld = 0; deltaw = 0; };


 protected:

//===========================================================================
/*!
 *  \brief Determines the sign of "x".
 *
 *  \param x The value of which the sign shall be determined.
 *  \return "-1", if the sign of \em x is negative, "0" otherwise.
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
  double    sgn(double x) 
  { 
    if (x > 0) return 1.; if (x < 0) return -1.; return 0.; 
  }

  //! The update value for the currently processed weight. 
  double delta;

  //! The update values for all weights.
  Array<double> deltaw;     

  //! The last error gradient.
  Array<double> dedwOld;

};

#endif

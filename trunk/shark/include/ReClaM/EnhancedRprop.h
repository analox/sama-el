//===========================================================================
/*!
 *  \file EnhancedRprop.h
 *
 *  \brief This file offers classes to use the 
 *         Resilient-Backpropagation-algorithm for the optimization of the
 *         adaptive parameters of a network.
 *         
 *  Four classes with four versions of the algorithm are included in this
 *  file: <br>
 *
 *  <ul>
 *      <li>IRpropMinusTrigPower: An improved Rprop algorithm without 
 *                       weight-backtracking enhanced with Trig and Power.</li>
 *      <li>IRpropMinusPower: An improved Rprop algorithm without 
 *                       weight-backtracking enhanced with Power.</li>
 *      
 *  </ul>
 *
 *  \author  B. Mersch, C. Igel
 *  \date    2004
 *
 *  \par Copyright (c) 2004-2000:
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
 *      $RCSfile: EnhancedRprop.h,v $<BR>
 *      $Revision: 2.2 $<BR>
 *      $Date: 2005/08/03 12:04:30 $
 *
 *  \par Changes:
 *      $Log: EnhancedRprop.h,v $
 *      Revision 2.2  2005/08/03 12:04:30  christian_igel
 *      changes limits/values, MAXDOUBLE etc.
 *
 *      Revision 2.1  2004/06/02 16:44:25  shark-admin
 *      documentation enhanced
 *
 *      Revision 2.0  2004/05/04 13:39:27  saviapbe
 *
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.2  2004/05/04 13:37:13  saviapbe
 *
 *       CVS logging was set correctly and the explicit copy constructor was deleted.
 *

 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
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

#ifndef ENHANCED_RPROP_H
#define ENHANCED_RPROP_H

#include <Array/Array.h>
#include <ReClaM/ModelInterface.h>

#ifdef __SOLARIS__
#include <values.h>
#else
#include <limits>
#endif

//===========================================================================
/*!
 *  \brief This class offers methods for the usage of the improved
 *         Resilient-Backpropagation-algorithm without weight-backtracking
 *         enhanced with Trig and Power.
 */

class IRpropMinusTrigPower : virtual public ModelInterface {
 public:

//===========================================================================
/*!
 *  \brief Prepares the Rprop algorithm for the currently used network.
 *
 *  Internal variables of the class instance are initialized and memory
 *  for used structures adapted to the used network. An initial value
 *  for the parameter \f$\Delta\f$ is assigned to all weights of the network.
 *
 *  \param _delta0 Initial value for the parameter \f$\Delta\f$. 
 *                 The default value is "0.01".
 *  \return none
 *
 *  \author  B. Mersch
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void initRprop(double _delta0 = 0.0000000001) 
  { 
    deltaw.resize(w.nelem());
    delta.resize (w.nelem());
    dedwOld.resize (w.nelem());
    delta  = _delta0; 
    deltaw = 0; 
    dedwOld  = 0; 
    delta0 = _delta0;
  }


//===========================================================================
/*!
 *  \brief Performs iteration of the Rprop algorithm.
 *
 *  The error \f$E^{(t)}\f$ of the used network for the current iteration
 *  \f$t\f$ is calculated and the values of the weights \f$w_i,\ i = 1, 
 *  \dots, n\f$ and the parameters \f$\Delta_i\f$ are 
 *  adapted depending on this error.
 *
 *  \param in    The input patterns used for the training of the network.
 *  \param out   The target values for to the input patterns.
 *  \param np    The increase factor \f$\eta^+\f$, by default set to "1.2".
 *  \param nm    The decrease factor \f$\eta^-\f$, by default set to "0.5".
 *  \param dMax  The upper limit of the increments \f$\Delta w_i^{(t)}\f$.
 *               The default value is "50".
 *  \param dMin  The lower limit of the increments \f$\Delta w_i^{(t)}\f$.
 *               The default value is "1e-6".
 *  \return none
 *
 *  \author  C. Igel
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
 void rprop(const Array<double> &in,
	     const Array<double> &out,
	     double              np   = 1.05  ,
	     double              nm   = 0  ,
	     double              dMax = 1   ,
	     double              dMin = 1e-6 ) 

{
    derror(in, out, false);

    for(unsigned  i = 0; i < w.nelem(); i++) {
     
      if(dedw(i) * dedwOld(i) > 0) {
	delta(i) = pow(delta(i),(1/np));
	deltaw(i) = dMax * delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) = dedw(i);
      } else if(dedw(i) * dedwOld(i) < 0) {
	double  nm = cos(atan(dedw(i)))/(cos(atan(dedwOld(i))) + cos(atan(dedw(i))));
	delta(i) = nm * delta(i);
	dedwOld(i) = 0;
      } else {
	deltaw(i) = dMax * delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) = dedw(i);
      }
    }
  }
  


//===========================================================================
/*!
 *  \brief Resets internal variables of the class instance.
 *
 *      \return 
 *          none
 *
 *  \author  B. Mersch, C. Igel
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void clearMemory() { delta = delta0; dedwOld = 0; deltaw = 0; };


 protected:

//===========================================================================
/*!
 *  \brief Determines the sign of "x".
 *
 *  \param x The value of which the sign shall be determined.
 *  \return "-1", if the sign of \em x is negative, "0" otherwise.
 *
 *  \author  C. Igel
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
  int    sgn(double x) 
  { 
    if (x > 0) return 1; if (x < 0) return -1; return 0; 
  }

  //! Used to replace the standard definition.
  double min(double x, double y) { return x < y ? x : y; }

  //! Used to replace the standard definition.
  double max(double x, double y) { return x > y ? x : y; }

  //! The initial increment for all weights.
  double delta0;

  //! The final update values for all weights.
  Array<double> deltaw;     

  //! The last error gradient.
  Array<double> dedwOld;

  //! The absolute update values (increment) for all weights.
  Array<double> delta;

};

//===========================================================================
/*!
 *  \brief This class offers methods for the usage of the improved
 *         Resilient-Backpropagation-algorithm without weight-backtracking
 *         enhanced with Power.
 */

class IRpropMinusPower : virtual public ModelInterface {
 public:


#ifndef EXPLICIT_COPY
  // some compilers require this
  IRpropMinusPower &operator=(const IRpropMinusPower &net) {
    deltaw  = net.deltaw;
    delta   = net.delta;
    dedwOld = net.dedwOld;

    delta0  = net.delta0;

    return *this;
  }
#endif

//===========================================================================
/*!
 *  \brief Prepares the Rprop algorithm for the currently used network.
 *
 *  Internal variables of the class instance are initialized and memory
 *  for used structures adapted to the used network. An initial value
 *  for the parameter \f$\Delta\f$ is assigned to all weights of the network.
 *
 *  \param _delta0 Initial value for the parameter \f$\Delta\f$. 
 *                 The default value is "0.01".
 *  \return none
 *
 *  \author  B. Mersch, C. Igel
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void initRprop(double _delta0 = 0.0000000001) 
  { 
    deltaw.resize(w.nelem());
    delta.resize (w.nelem());
    dedwOld.resize (w.nelem());
    delta  = _delta0; 
    deltaw = 0; 
    dedwOld  = 0; 
    delta0 = _delta0;
  }


//===========================================================================
/*!
 *  \brief Performs iteration of the Rprop algorithm.
 *
 *  The error \f$E^{(t)}\f$ of the used network for the current iteration
 *  \f$t\f$ is calculated and the values of the weights \f$w_i,\ i = 1, 
 *  \dots, n\f$ and the parameters \f$\Delta_i\f$ are 
 *  adapted depending on this error.
 *
 *  \author  B. Mersch, C. Igel
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void rprop(const Array<double> &in,
	     const Array<double> &out,
	     double              np   = 1.05  ,
	     double              nm   = .5   ,
	     double              dMax = 1   ,
	     double              dMin = 1e-6 ) 
  
{
    derror(in, out, false);

    for(unsigned  i = 0; i < w.nelem(); i++) {
      if(dedw(i) * dedwOld(i) > 0) {
	delta(i) = pow(delta(i),(1/np)); 
	deltaw(i) = dMax * delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) = dedw(i);
      } else if(dedw(i) * dedwOld(i) < 0) {
	delta(i) = nm * delta(i);
	dedwOld(i) = 0;
      } else {
	deltaw(i) = dMax * delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) =  dedw(i);
      }
    }
  }


//===========================================================================
/*!
 *  \brief Resets internal variables of the class instance.
 *
 *      \return 
 *          none
 *
 *  \author  B. Mersch, C. Igel
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void clearMemory() { delta = delta0; dedwOld = 0; deltaw = 0; };


 protected:

//===========================================================================
/*!
 *  \brief Determines the sign of "x".
 *
 *  \param x The value of which the sign shall be determined.
 *  \return "-1", if the sign of \em x is negative, "0" otherwise.
 *
 *  \author  C. Igel
 *  \date    2004
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
  int    sgn(double x) 
  { 
    if (x > 0) return 1; if (x < 0) return -1; return 0; 
  }

  //! Used to replace the standard definition.
  double min(double x, double y) { return x < y ? x : y; }

  //! Used to replace the standard definition.
  double max(double x, double y) { return x > y ? x : y; }

  //! The initial increment for all weights.
  double delta0;

  //! The final update values for all weights.
  Array<double> deltaw;     

  //! The last error gradient.
  Array<double> dedwOld;

  //! The absolute update values (increment) for all weights.
  Array<double> delta;

};
#endif

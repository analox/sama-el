//===========================================================================
/*!
 *  \file Rprop.h
 *
 *  \brief This file offers classes to use the 
 *         Resilient-Backpropagation-algorithm for the optimization of the
 *         adaptive parameters of a network.
 *         
 *  Four classes with four versions of the algorithm are included in this
 *  file: <br>
 *
 *  <ul>
 *      <li>RpropPlus:   The Rprop algorithm with weight-backtracking.</li>
 *      <li>RpropMinus:  The Rprop algorithm without weight-backtracking.</li>
 *      <li>IRpropPlus:  An improved Rprop algorithm with 
 *                       weight-backtracking.</li>
 *      <li>IRpropMinus: An improved Rprop algorithm without 
 *                       weight-backtracking.</li>
 *  </ul>
 *
 *  \author  C. Igel, M. H&uuml;sken
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
 *      $RCSfile: Rprop.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2005/08/03 12:02:33 $
 *
 *  \par Changes:
 *      $Log: Rprop.h,v $
 *      Revision 2.1  2005/08/03 12:02:33  christian_igel
 *      changed limits/values etc.
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
 *      Revision 1.3  2002/02/06 14:48:29  rudi
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

#ifndef RPROP_H
#define RPROP_H

#include "Array/Array.h"
#include "ReClaM/ModelInterface.h"
#ifndef __SOLARIS__
#include <limits>
#else
#include <values.h>
#endif


//===========================================================================
/*!
 *  \brief This class offers methods for the usage of the 
 *         Resilient-Backpropagation-algorithm with weight-backtracking.
 *
 *  The Rprop algorithm is an improvement of the algorithms with adaptive
 *  learning rates (as the Adaptive Backpropagation algorithm by Silva
 *  and Ameida, please see AdpBP.h for a description of the 
 *  working of such an algorithm), that uses increments for the update
 *  of the weights, that are independant from the absolute partial
 *  derivatives. This makes sense, because large flat regions
 *  in the search space (plateaus) lead to small absolute partial
 *  derivatives and so the increments are chosen small, but the increments
 *  should be large to skip the plateau. In contrast, the absolute partial
 *  derivatives are very large at the "slopes" of very "narrow canyons",
 *  which leads to large increments that will skip the minimum lying
 *  at the bottom of the canyon, but it would make more sense to
 *  chose small increments to hit the minimum. <br>
 *  So, the Rprop algorithm only uses the signs of the partial derivatives
 *  and not the absolute values to adapt the parameters. <br>
 *  Instead of individual learning rates, it uses the parameter
 *  \f$\Delta_i^{(t)}\f$ for weight \f$w_i,\ i = 1, \dots, n\f$ in 
 *  iteration \f$t\f$, where the parameter will be adapted before the 
 *  change of the weights: <br>
 *
 *  \f$
 *  \Delta_i^{(t)} = \Bigg\{ 
 *  \begin{array}{ll}
 *  min( \eta^+ \cdot \Delta_i^{(t-1)}, \Delta_{max} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} > 0 \\
 *  max( \eta^- \cdot \Delta_i^{(t-1)}, \Delta_{min} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0 \\
 *  \Delta_i^{(t-1)}, & \mbox{otherwise}
 *  \end{array}
 *  \f$ 
 *
 *  The parameters \f$\eta^+ > 1\f$ and \f$0 < \eta^- < 1\f$ control
 *  the speed of the adaptation. To stabilize the increments, they are
 *  restricted to the interval \f$[\Delta_{min}, \Delta_{max}]\f$. <br>
 *  After the adaptation of the \f$\Delta_i\f$ the update for the 
 *  weights will be calculated as
 *
 *  \f$
 *  \Delta w_i^{(t)} := - \mbox{sign} 
 *  \left( \frac{\partial E^{(t)}}{\partial w_i}\right) \cdot \Delta_i^{(t)}
 *  \f$
 *        
 *  Furthermore, weight-backtracking will take place to increase the
 *  stability of the method, i.e.
 *  if \f$\frac{\partial E^{(t-1)}}{\partial w_i} \cdot 
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0\f$ then
 *  \f$\Delta w_i^{(t)} := - \Delta w_i^{(t-1)}; 
 *  \frac{\partial E^{(t)}}{\partial w_i} := 0\f$, where
 *  the assignment of zero to the partial derivative of the error
 *  leads to a freezing of the increment in the next iteration.
 *
 *  For further information about the algorithm, please refer to: <br>
 *  
 *  Martin Riedmiller and Heinrich Braun, <br>
 *  "A Direct Adaptive Method for Faster Backpropagation Learning: The 
 *  RPROP Algorithm". <br> 
 *  In "Proceedings of the IEEE International Conference on Neural Networks",
 *  pp. 586-591, <br>
 *  Published by IEEE Press in 1993 
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
class RpropPlus : virtual public ModelInterface {
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
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void initRprop(double _delta0 = 0.01) 
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
 *  \brief Performs a run of the Rprop algorithm.
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
 *  \date    1999
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
	     double              np   = 1.2  ,
	     double              nm   = .5   ,
	     double              dMax = 50   ,
	     double              dMin = 1e-6 ) 
  {
    derror(in, out, false);

    for(unsigned  i = 0; i < w.nelem(); i++) {
      if(dedw(i) * dedwOld(i) > 0) {
	delta(i) = min(dMax, np * delta(i));
	deltaw(i) = delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) = dedw(i);
      } else if(dedw(i) * dedwOld(i) < 0) {
	delta(i) = max(dMin, nm * delta(i));
	w(i) -= deltaw(i);
	dedwOld(i) = 0;
      } else {
	deltaw(i) = delta(i) * -sgn(dedw(i));
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
 *  \author  C. Igel
 *  \date    1999
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
 *  \date    1999
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
 *  \brief This class offers methods for the usage of the 
 *         Resilient-Backpropagation-algorithm without weight-backtracking.
 *
 *  The Rprop algorithm is an improvement of the algorithms with adaptive
 *  learning rates (as the Adaptive Backpropagation algorithm by Silva
 *  and Ameida, please see AdpBP.h for a description of the 
 *  working of such an algorithm), that uses increments for the update
 *  of the weights, that are independant from the absolute partial
 *  derivatives. This makes sense, because large flat regions
 *  in the search space (plateaus) lead to small absolute partial
 *  derivatives and so the increments are chosen small, but the increments
 *  should be large to skip the plateau. In contrast, the absolute partial
 *  derivatives are very large at the "slopes" of very "narrow canyons",
 *  which leads to large increments that will skip the minimum lying
 *  at the bottom of the canyon, but it would make more sense to
 *  chose small increments to hit the minimum. <br>
 *  So, the Rprop algorithm only uses the signs of the partial derivatives
 *  and not the absolute values to adapt the parameters. <br>
 *  Instead of individual learning rates, it uses the parameter
 *  \f$\Delta_i^{(t)}\f$ for weight \f$w_i,\ i = 1, \dots, n\f$ in 
 *  iteration \f$t\f$, where the parameter will be adapted before the 
 *  change of the weights: <br>
 *
 *  \f$
 *  \Delta_i^{(t)} = \Bigg\{ 
 *  \begin{array}{ll}
 *  min( \eta^+ \cdot \Delta_i^{(t-1)}, \Delta_{max} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} > 0 \\
 *  max( \eta^- \cdot \Delta_i^{(t-1)}, \Delta_{min} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0 \\
 *  \Delta_i^{(t-1)}, & \mbox{otherwise}
 *  \end{array}
 *  \f$ 
 *
 *  The parameters \f$\eta^+ > 1\f$ and \f$0 < \eta^- < 1\f$ control
 *  the speed of the adaptation. To stabilize the increments, they are
 *  restricted to the interval \f$[\Delta_{min}, \Delta_{max}]\f$. <br>
 *  After the adaptation of the \f$\Delta_i\f$ the update for the 
 *  weights will be calculated as
 *
 *  \f$
 *  \Delta w_i^{(t)} := - \mbox{sign} 
 *  \left( \frac{\partial E^{(t)}}{\partial w_i}\right) \cdot \Delta_i^{(t)}
 *  \f$
 *        
 *  For further information about the algorithm, please refer to: <br>
 *  
 *  Martin Riedmiller, <br>
 *  "Advanced Supervised Learning in Multi-layer Perceptrons - 
 *  From Backpropagation to Adaptive Learning Algorithms". <br>
 *  In "International Journal of Computer Standards and Interfaces", volume 16,
 *  no. 5, 1994, pp. 265-278 <br> 
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
class RpropMinus : virtual public ModelInterface {
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
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void initRprop(double _delta0 = 0.01) 
  { 
    delta.resize (w.nelem());
    dedwOld.resize (w.nelem());
    delta  = _delta0; 
    dedwOld  = 0; 
    delta0 = _delta0;
  }


//===========================================================================
/*!
 *  \brief Performs a run of the Rprop algorithm.
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
 *  \date    1999
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
	     double              np   = 1.2  ,
	     double              nm   = .5   ,
	     double              dMax = 50   ,
	     double              dMin = 1e-6 )
  {
    derror(in, out,false);
    for(unsigned  i = 0; i < w.nelem(); i++) {
      if(dedw(i) * dedwOld(i) > 0) {
	delta(i) = min(dMax, np * delta(i));
      } else if(dedw(i) * dedwOld(i) < 0) {
	delta(i) = max(dMin, nm * delta(i));
      } else {
	;
      }
      w(i) += delta(i) * -sgn(dedw(i));
      dedwOld(i) = dedw(i);
    }
  }


//===========================================================================
/*!
 *  \brief Resets internal variables of the class instance.
 *
 *      \return 
 *          none
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
 *
 */  
  void clearMemory() { delta = delta0; dedwOld = 0; };


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
 *         Resilient-Backpropagation-algorithm with weight-backtracking.
 *
 *  The Rprop algorithm is an improvement of the algorithms with adaptive
 *  learning rates (as the Adaptive Backpropagation algorithm by Silva
 *  and Ameida, please see AdpBP.h for a description of the 
 *  working of such an algorithm), that uses increments for the update
 *  of the weights, that are independant from the absolute partial
 *  derivatives. This makes sense, because large flat regions
 *  in the search space (plateaus) lead to small absolute partial
 *  derivatives and so the increments are chosen small, but the increments
 *  should be large to skip the plateau. In contrast, the absolute partial
 *  derivatives are very large at the "slopes" of very "narrow canyons",
 *  which leads to large increments that will skip the minimum lying
 *  at the bottom of the canyon, but it would make more sense to
 *  chose small increments to hit the minimum. <br>
 *  So, the Rprop algorithm only uses the signs of the partial derivatives
 *  and not the absolute values to adapt the parameters. <br>
 *  Instead of individual learning rates, it uses the parameter
 *  \f$\Delta_i^{(t)}\f$ for weight \f$w_i,\ i = 1, \dots, n\f$ in 
 *  iteration \f$t\f$, where the parameter will be adapted before the 
 *  change of the weights: <br>
 *
 *  \f$
 *  \Delta_i^{(t)} = \Bigg\{ 
 *  \begin{array}{ll}
 *  min( \eta^+ \cdot \Delta_i^{(t-1)}, \Delta_{max} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} > 0 \\
 *  max( \eta^- \cdot \Delta_i^{(t-1)}, \Delta_{min} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0 \\
 *  \Delta_i^{(t-1)}, & \mbox{otherwise}
 *  \end{array}
 *  \f$ 
 *
 *  The parameters \f$\eta^+ > 1\f$ and \f$0 < \eta^- < 1\f$ control
 *  the speed of the adaptation. To stabilize the increments, they are
 *  restricted to the interval \f$[\Delta_{min}, \Delta_{max}]\f$. <br>
 *  After the adaptation of the \f$\Delta_i\f$ the update for the 
 *  weights will be calculated as
 *
 *  \f$
 *  \Delta w_i^{(t)} := - \mbox{sign} 
 *  \left( \frac{\partial E^{(t)}}{\partial w_i}\right) \cdot \Delta_i^{(t)}
 *  \f$
 *        
 *  Furthermore, weight-backtracking will take place to increase the
 *  stability of the method. In contrast to the original Rprop algorithm
 *  with weight-backtracking (see RpropPlus) this weight-backtracking 
 *  is improved by additionally taken the error of the last iteration
 *  \f$t - 1\f$ into account. <br>
 *  The idea of this modification is, that a change of the sign of the
 *  partial derivation \f$\frac{\partial E}{\partial w_i}\f$
 *  only states, that a minimum was skipped and not, whether this step
 *  lead to an approach to the minimum or not. <br>
 *  By using the old error value the improved weight-backtracking only
 *  undoes changes, when the error has increased and only the parameters
 *  \f$w_i\f$ are reset to the old values, where a sign change of
 *  \f$\frac{\partial E}{\partial w_i}\f$ has taken place. <br>
 *  So the new weight-backtracking rule is: <br>
 *
 *  \f$
 *  \mbox{if\ } \frac{\partial E^{(t-1)}}{\partial w_i} \cdot 
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0 \mbox{\ then} \{
 *  \f$
 *
 *  \f$  
 *  \begin{array}{lll}
 *   \Delta w_i^{(t)} = \bigg\{ &
 *   - \Delta w_i^{(t-1)}, & \mbox{if\ } E^{(t)} > E^{(t - 1)} \\
 *   & 0, & otherwise \\
 *  \frac{\partial E^{(t)}}{\partial w_i} := 0 
 *  \end{array}
 *  \f$
 *
 *  \f$\}\f$
 *
 *  , where the assignment of zero to the partial derivative of the error
 *  leads to a freezing of the increment in the next iteration. <br>
 *
 *  This modification of the weight backtracking leads to a better
 *  optimization on artifical, paraboloidal error surfaces. <br>
 *
 *  For further information about the algorithm, please refer to: <br>
 *  
 *  Christian Igel and Michael H&uuml;sken, <br>
 *  "Empirical Evaluation of the Improved Rprop Learning Algorithm". <br>
 *  In Neurocomputing Journal, 2002, in press <br> 
 *
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 *
 */
class IRpropPlus : virtual public ModelInterface {
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
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void initRprop(double _delta0 = 0.01) 
  { 
    deltaw.resize(w.nelem());
    delta.resize (w.nelem());
    dedwOld.resize (w.nelem());
    delta  = _delta0; 
    deltaw = 0; 
    dedwOld  = 0; 
    delta0 = _delta0;
#ifndef __SOLARIS__
    oldE   = std::numeric_limits< double >::max( );
#else
    oldE   = MAXDOUBLE;
#endif
  }


//===========================================================================
/*!
 *  \brief Performs a run of the Rprop algorithm.
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
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  double rprop(const Array<double> &in,
	     const Array<double> &out,
	     double              np   = 1.2  ,
	     double              nm   = .5   ,
	     double              dMax = 50   ,
	     double              dMin = 1e-6 )
  {
    double newE = derror(in, out);

    for(unsigned  i = 0; i < w.nelem(); i++) {
      if(dedw(i) * dedwOld(i) > 0) {
	delta(i) = min(dMax, np * delta(i));
	deltaw(i) = delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) = dedw(i);
      } else if(dedw(i) * dedwOld(i) < 0) {
	delta(i) = max(dMin, nm * delta(i));
	//cout << "Oho!" << endl;
	if(oldE < newE) {
	  //cout << "Aha! ----------------" << endl;
	  w(i) -= deltaw(i);
	}
	dedwOld(i) = 0;
      } else {
	deltaw(i) = delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) = dedw(i);
      }
    }
    oldE = newE;
    return newE; 
  }


//===========================================================================
/*!
 *  \brief Resets internal variables of the class instance.
 *
 *      \return 
 *          none
 *
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void clearMemory()
  {
	  delta = delta0;
	  dedwOld = 0;
	  deltaw = 0;
#ifndef __SOLARIS__
	  oldE = std::numeric_limits< double >::max( );
#else
	  oldE = MAXDOUBLE;
#endif
  }


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

  //! The error of the last iteration.
  double oldE;

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
 *         Resilient-Backpropagation-algorithm without weight-backtracking.
 *
 *  The Rprop algorithm is an improvement of the algorithms with adaptive
 *  learning rates (as the Adaptive Backpropagation algorithm by Silva
 *  and Ameida, please see AdpBP.h for a description of the 
 *  working of such an algorithm), that uses increments for the update
 *  of the weights, that are independant from the absolute partial
 *  derivatives. This makes sense, because large flat regions
 *  in the search space (plateaus) lead to small absolute partial
 *  derivatives and so the increments are chosen small, but the increments
 *  should be large to skip the plateau. In contrast, the absolute partial
 *  derivatives are very large at the "slopes" of very "narrow canyons",
 *  which leads to large increments that will skip the minimum lying
 *  at the bottom of the canyon, but it would make more sense to
 *  chose small increments to hit the minimum. <br>
 *  So, the Rprop algorithm only uses the signs of the partial derivatives
 *  and not the absolute values to adapt the parameters. <br>
 *  Instead of individual learning rates, it uses the parameter
 *  \f$\Delta_i^{(t)}\f$ for weight \f$w_i,\ i = 1, \dots, n\f$ in 
 *  iteration \f$t\f$, where the parameter will be adapted before the 
 *  change of the weights. <br>
 *  As an improving modification, this algorithm
 *  adapts the "freezing" of the increment in the next iteration as
 *  usually only practiced by the Rprop algorithm with weight-backtracking
 *  (see RpropPlus), i.e. \f$\frac{\partial E^{(t)}}{\partial w_i} := 0\f$.
 *  Tests have shown a far more better optimization when using this 
 *  modification. So the new adaptation rule of \f$\Delta\f$ is given
 *  as: <br>
 *
 *  \f$
 *  \Delta_i^{(t)} = \Bigg\{ 
 *  \begin{array}{ll}
 *  min( \eta^+ \cdot \Delta_i^{(t-1)}, \Delta_{max} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} > 0 \\
 *  max( \eta^- \cdot \Delta_i^{(t-1)}, \Delta_{min} ); 
 *  \frac{\partial E^{(t)}}{\partial w_i} := 0, & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0 \\
 *  \Delta_i^{(t-1)}, & \mbox{otherwise}
 *  \end{array}
 *  \f$ 
 *
 *  The parameters \f$\eta^+ > 1\f$ and \f$0 < \eta^- < 1\f$ control
 *  the speed of the adaptation. To stabilize the increments, they are
 *  restricted to the interval \f$[\Delta_{min}, \Delta_{max}]\f$. <br>
 *  After the adaptation of the \f$\Delta_i\f$ the update for the 
 *  weights will be calculated as
 *
 *  \f$
 *  \Delta w_i^{(t)} := - \mbox{sign} 
 *  \left( \frac{\partial E^{(t)}}{\partial w_i}\right) \cdot \Delta_i^{(t)}
 *  \f$
 *        
 *  For further information about the algorithm, please refer to: <br>
 *
 *  Christian Igel and Michael H&uuml;sken, <br>
 *  "Empirical Evaluation of the Improved Rprop Learning Algorithm". <br>
 *  In Neurocomputing Journal, 2002, in press <br> 
 *
 *
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 *
 */
class IRpropMinus : virtual public ModelInterface {
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
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void initRprop(double _delta0 = 0.01) 
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
 *  \brief Performs a run of the Rprop algorithm.
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
 *  \date    1999
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
	     double              np   = 1.2  ,
	     double              nm   = .5   ,
	     double              dMax = 50   ,
	     double              dMin = 1e-6 ) 
  {
    derror(in, out, false);

    for(unsigned  i = 0; i < w.nelem(); i++) {
      if(dedw(i) * dedwOld(i) > 0) {
	delta(i) = min(dMax, np * delta(i));
	deltaw(i) = delta(i) * -sgn(dedw(i));
	w(i) += deltaw(i);
	dedwOld(i) = dedw(i);
      } else if(dedw(i) * dedwOld(i) < 0) {
	delta(i) = max(dMin, nm * delta(i));
	dedwOld(i) = 0;
      } else {
	deltaw(i) = delta(i) * -sgn(dedw(i));
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
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
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
 *  \date    1999
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

//! Used to connect the class names with the year of 
//! publication of the paper in which the algorithm was introduced.
typedef IRpropPlus Rprop99;

//! Used to connect the class names with the year of 
//! publication of the paper in which the algorithm was introduced.
typedef IRpropMinus Rprop99d;

//! Used to connect the class names with the year of 
//! publication of the paper in which the algorithm was introduced.
typedef RpropPlus Rprop93;

//! Used to connect the class names with the year of 
//! publication of the paper in which the algorithm was introduced.
typedef RpropMinus Rprop94;

#endif





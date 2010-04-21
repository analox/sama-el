//===========================================================================
/*!
 *  \file AdpBP.h
 *
 *  \brief Offers the two versions of the gradient descent based
 *         optimization algorithm with individual adaptive learning
 *         rates by Silva and Almeida (Adaptive BackPropagation).
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
 *      C++
 *
 *  \par File and Revision:
 *      $RCSfile: AdpBP.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2005/08/03 12:02:33 $
 *
 *  \par Changes:
 *      $Log: AdpBP.h,v $
 *      Revision 2.1  2005/08/03 12:02:33  christian_igel
 *      changed limits/values etc.
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:24:13  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2002/02/06 14:30:12  rudi
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


#ifndef ADP_BP_H
#define ADP_BP_H

#include "Array/ArrayOp.h"
#include "ReClaM/ModelInterface.h"
#ifndef __SOLARIS__
#include <limits>
#else
#include <values.h>
#endif


//===========================================================================
/*!
 *  \brief Offers the gradient-based optimization algorithm with 
 *         individual adaptive learning rates by Silva and Almeida
 *         (Adaptive BackPropagation).
 *
 *  This optimization algorithm introduced by Silva and Almeida adds
 *  individual adaptive learning rates and weight-backtracking to
 *  standard steepest descent.  
 *  \par 
 *  To avoid bad choices of the learning rate, the algorithm by
 *  Silva and Almeida uses an individual adaptive learning rate
 *  \f$\eta_i\f$ for each weight \f$w_i\f$.  The adaptation of
 *  \f$\eta_i\f$ is determined by the sign of
 *  the partial derivation \f$\frac{\partial E}{\partial w_i}\f$. If
 *  the sign of the derivation changes from iteration \f$t - 1\f$ to
 *  iteration \f$t\f$, then at least one minimum with regard to
 *  \f$w_i\f$ was skipped and hence the learning rate \f$\eta_i\f$ is
 *  decreased by the multiplication with a predefined constant.
 *  If there is no change of the sign, it a plateau of the objective function
 *  is assumed and the learning rate \f$\eta_i\f$ is increased by the
 *  multiplication with a second constant.  Furthermore,
 *  weight-backtracking is used, i.e.. if the error increases from
 *  iteration \f$t - 1\f$ to \f$t\f$, then the last weight adaptation
 *  is undone. For further information about this algorithm,
 *  please refer to: 
 *
 \verbatim
 @InCollection{silva:90,
    author =       {Fernando M. Silva and Luis B. Almeida},
    editor =       {Rolf Eckmiller},
    booktitle =    {Advanced Neural Computers},
    title =        {Speeding Up Backpropagation},
    publisher =    {North-Holland},
    year =         {1990},
    pages =        {151-158},
  }
 \endverbatim
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable, documentation checked
 *
 *  \sa AdpBP90b
 *
 */
class AdpBP90a : virtual public ModelInterface {

 public:


//===========================================================================
/*!
 *  \brief Prepares the optimization algorithm for the currently used network.
 *
 *  Internal variables of the class instance are initialized. An initial learning
 *  rate \f$\eta_0\f$ is assigned to all weights of the network.
 *
 *  \param _eta0 Initial value for the adaptive learning rate \f$\eta\f$
 *               for all values. The default value is "0.1".
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
  void initAdpBP(double _eta0 = 0.1) 
  { 
    deltaw.resize(w.nelem());
    eta   .resize(w.nelem());
    dedwOld .resize(w.nelem());
    v     .resize(w.nelem());
    v      = 0.;
    eta    = _eta0; 
    deltaw = 0; 
    dedwOld  = 0; 
    eta0   = _eta0;
#ifndef __SOLARIS__
    oldE   = std::numeric_limits< double >::max( );
#else
    oldE   = MAXDOUBLE;
#endif
  }



//===========================================================================
/*!
 *  \brief Performs a run of the optimization algorithm.
 *
 *  The error \f$E^{(t)}\f$ and its derivatives with respect to the model 
 *  parameters for the current iteration 
 *  \f$t\f$ are calculated and the values of the weights \f$w_i,\ i = 1, 
 *  \dots, n\f$ and the  individual learning rates \f$\eta_i\f$ are 
 *  adapted.
 *
 *  \param in    The input patterns used for the training of the model.
 *  \param out   The target values for the input patterns.
 *  \param u     The increase factor for the adaptive learning rate.
 *               The default value is "1.2".
 *  \param d     The decrease factor for the adaptive learning rate.
 *               The default value is "0.7".
 *  \param alpha The momentum factor \f$\alpha\f$ that controls the effect
 *               of the momentum term, by default set to "0.5".
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
  void adpBP(const Array<double> &in,
	     const Array<double> &out,
	     double              u     = 1.2,
	     double              d     = .7,
	     double              alpha = 0.5)   
  {
    double newE = derror(in, out);

    if(oldE < newE) w -= deltaw;

    for(unsigned  i = 0; i < w.nelem(); i++) {
      if(dedw(i) * dedwOld(i) > 0) eta(i) *= u;
      if(dedw(i) * dedwOld(i) < 0) eta(i) *= d;
    }

    v = dedw + alpha * v;
    deltaw = -eta * v; 

    w += deltaw;

    dedwOld = dedw;
    oldE = newE;
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
  void clearMemory() 
  { 
#ifndef __SOLARIS__
    oldE   = std::numeric_limits< double >::max( );
#else
    oldE   = MAXDOUBLE;
#endif
    eta = eta0; dedwOld = 0; deltaw = 0;
  };


 protected:

  //! The initial value for the adaptive learning rate \f$\eta\f$
  //! for all weights of the network.
  double eta0;

  //! The old error value.
  double oldE;

  //! The final update values (including the individual learning rates)
  //! for all weights.
  Array<double> deltaw;

  //! The update values (without the individual learning rates) for all
  //! weights.
  Array<double> v;     

  //! The last error gradient.
  Array<double> dedwOld;

  //! The adaptive learning rate \f$\eta\f$
  //! (for each weight individually).
  Array<double> eta;

};


//===========================================================================
/*!
 *  \brief Offers the second version of the gradient descent based 
 *         optimization algorithm with individual adaptive learning rates 
 *         by Silva and Almeida (Adaptive BackPropagation).
 *  This optimization algorithm introduced by Silva and Almeida adds
 *  individual adaptive learning rates and weight-backtracking to
 *  standard steepest descent.  
 *  \par 
 *  To avoid bad choices of the learning rate, the algorithm by
 *  Silva and Almeida uses an individual adaptive learning rate
 *  \f$\eta_i\f$ for each weight \f$w_i\f$.  The adaptation of
 *  \f$\eta_i\f$ is determined by the sign of
 *  the partial derivation \f$\frac{\partial E}{\partial w_i}\f$. If
 *  the sign of the derivation changes from iteration \f$t - 1\f$ to
 *  iteration \f$t\f$, then at least one minimum with regard to
 *  \f$w_i\f$ was skipped and hence the learning rate \f$\eta_i\f$ is
 *  decreased by the multiplication with a predefined constant.
 *  If there is no change of the sign, it a plateau of the objective function
 *  is assumed and the learning rate \f$\eta_i\f$ is increased by the
 *  multiplication with a second constant.  Furthermore,
 *  weight-backtracking is used, i.e.. if the error increases from
 *  iteration \f$t - 1\f$ to \f$t\f$, then the last weight adaptation
 *  is undone. For further information about this algorithm,
 *  please refer to: 
 *
 \verbatim
 @InCollection{silva:90b,
  author =       {Fernando M. Silva and Luis B. Almeida},
  title =        {Acceleration Techniques for the Backpropagation Algorithm},
  booktitle =    {Neural Networks -- EURASIP Workshop 1990},
  pages =        {110-119},
  year =         {1990},
  editor =       {Luis B. Almeida and C. J. Wellekens},
  number =       {412},
  series =       {LNCS},
  publisher = {Springer-Verlag},
 }
 \endverbatim 
 *  This second version of the algorithm has some minor modifications
 *  that improve the performance. 
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
class AdpBP90b : virtual public ModelInterface {
 public:


//===========================================================================
/*!
 *  \brief Prepares the optimization algorithm for the currently used model.
 *
 *  Internal variables of the class instance are initialized. An
 *  initial learning rate \f$\eta_0\f$ is assigned to all weights of
 *  the network.
 *
 *  \param _eta0 Initial value for the adaptive learning rate \f$\eta\f$
 *               for all values. The default value is "0.1".
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
  void initAdpBP(double _eta0 = 0.1) 
  { 
    deltaw.resize(w.nelem());
    eta   .resize(w.nelem());
    dedwOld .resize(w.nelem());
    v     .resize(w.nelem());
    v      = 0.;
    eta    = _eta0; 
    deltaw = 0; 
    dedwOld  = 0; 
    eta0   = _eta0;
#ifndef __SOLARIS__
    oldE   = std::numeric_limits< double >::max( );
#else
    oldE   = MAXDOUBLE;
#endif
    eta = eta0; dedwOld = 0; deltaw = 0;
  };
}


//===========================================================================
/*!
 *  \brief Performs a run of the optimization algorithm.
 *
 *  The error \f$E^{(t)}\f$ and its derivatives with respect to the model 
 *  parameters for the current iteration 
 *  \f$t\f$ are calculated and the values of the weights \f$w_i,\ i = 1, 
 *  \dots, n\f$ and the  individual learning rates \f$\eta_i\f$ are 
 *  adapted.
 *
 *  \param in    The input patterns used for the training of the network.
 *  \param out   The target values for the input patterns.
 *  \param u     The increase factor for the adaptive learning rate.
 *               The default value is "1.2".
 *  \param d     The decrease factor for the adaptive learning rate.
 *               The default value is "0.7".
 *  \param alpha The momentum factor \f$\alpha\f$ that controls the effect
 *               of the momentum term, by default set to "0.5".
 *  \return none
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      2002-01-18, ra: <br>
 *      Removed some forgotten "couts" used for monitoring.
 *
 *  \par Status
 *      stable
 *
 *
 */  
  void adpBP(const Array<double> &in,
	     const Array<double> &out,
	     double              u     = 1.2  ,   // increase factor
	     double              d     = .7   ,   // decrease factor
	     double              alpha = 0.5) {   // momentum factor
    double newE = derror(in, out);


    for(unsigned  i = 0; i < w.nelem(); i++) {
      if(dedw(i) * dedwOld(i) > 0) eta(i) *= u;
      if(dedw(i) * dedwOld(i) < 0) eta(i) *= d;
    }

    if((1E-4 + oldE) < newE) {
      w -= deltaw;
      //deltaw = 0.;
      dedwOld  = dedw;
    }
    else {
      v      = dedw + alpha * v;
      deltaw = -eta * v; 
      dedwOld  = dedw;
      oldE   = newE;
      w += deltaw;
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
  void clearMemory() { 
    eta = eta0; dedwOld = 0; deltaw = 0;
#ifndef __SOLARIS__
    oldE   = std::numeric_limits< double >::max( );
#else
    oldE   = MAXDOUBLE;
#endif
  };


 protected:

  //! The initial value for the adaptive learning rate \f$\eta\f$
  //! for all weights of the network.
  double eta0;

  //! The old error value.
  double oldE;

  //! The final update values (including the individual learning rates)
  //! for all weights.
  Array<double> deltaw;

  //! The update values (without the individual learning rates) for all
  //! weights.
  Array<double> v;     

  //! The last error gradient.
  Array<double> dedwOld;

  //! The adaptive learning rate \f$\eta\f$
  //! (for each weight individually).
  Array<double> eta;

};

#endif



//===========================================================================
/*!
 *  \file CG.h
 *
 *  \brief Offers the Conjugate Gradients algorithm
 *         for the optimization of models. 
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Copyright (c) 1999
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
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
 *      $RCSfile: CG.h,v $<BR>
 *      $Revision: 2.2 $<BR>
 *      $Date: 2006/02/02 12:01:40 $
 *
 *  \par Changes:
 *      $Log: CG.h,v $
 *      Revision 2.2  2006/02/02 12:01:40  christian_igel
 *      little errors corrected, namespace changed
 *
 *      Revision 2.1  2005/08/03 12:07:48  christian_igel
 *      values/limits etc.
 *
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/02/06 14:30:46  rudi
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

#ifndef CG_H
#define CG_H

#include "ReClaM/ModelInterface.h"
#ifdef __SOLARIS__
	#include <values.h>
#else
	#include <limits>
#endif


//===========================================================================
/*!
 *  \brief Offers methods to use the Conjugate Gradients
 *         algorithm for the optimization of models.
 *
 *  As usual with other quasi-Newton methods, the error function
 *  is optimized by \f$n\f$ line searches, but here a special
 *  method is used to avoid bad performance that can occurr 
 *  when using simple successive gradient vectors. <br>
 *  When performing a line search from point \f$w^{(t)}\f$ along
 *  direction \f$d^{(t)}\f$, where the error minimum along
 *  this search path lies at \f$w^{(t + 1)}\f$ it holds (using
 *  the line search definitions): <br>
 *
 *  \f$
 *      w^{(t + 1)} = w^{(t)} + \lambda^{(t)} d^{(t)}, \mbox{\ with\ }
 *      \lambda \mbox{\ used for minimization, i.e.}
 *  \f$
 *
 *  \f$
 *      E_{min} = E(\lambda) = E( w^{(t)} + \lambda d^{(t)}).
 *  \f$
 *
 *  Looking at the gradient \f$g\f$ it follows:
 * 
 *  \f$
 *      \frac{\partial}{\partial \lambda} E(w^{(t)} + \lambda d^{(t)}) = 0
 *  \f$
 *
 *  which gives
 *
 *  \f$
 *      g^{(t + 1)T} d^{(t)} = 0, \mbox{\ where\ } g \equiv \nabla E.
 *  \f$
 *
 *  Thus, the gradient at the new minimum is orthogonal to the
 *  previous search direction. If one now chooses successive search
 *  directions to be the local (negative) gradient directions,
 *  this can lead to the search point oscillating on these
 *  successive steps while making only little progress towards the
 *  minimum - the minimization proceeds very slowly. <br> 
 *  A solution is to choose the next search direction \f$d^{(t + 1)}\f$
 *  such that, along this new direction, we retain the property
 *  that the component of the gradient parallel to the previous
 *  search direction remains zero (to lowest order): <br>
 *
 *  \f$
 *      g(w^{(t + 1)} + \lambda d^{(t + 1)})^Td^{(t)} = 0
 *  \f$
 *  
 *  When expanding to first order in \f$\lambda\f$, we finally 
 *  obtain: <br>
 *
 *  \f$
 *    d^{(t + 1)T} H d^{(t)} = 0
 *  \f$
 *
 *  where \f$H\f$ is the Hessian matrix evaluated at point \f$w^{(t + 1)}\f$. 
 *  <br>
 *  Search directions which satisfy this equation are said to be conjugate.
 *  For the evaluation of the search directions for quadratic target
 *  functions, the Polak-Ribi�e formula is used: <br>
 *
 *  \f$
 *    d^{( t + 1)} := - \nabla E \arrowvert_{w^{(t + 1)}} +
 *    \frac{\nabla E \arrowvert^T_{w^{(t + 1)}} (\nabla E 
 *    \arrowvert_{w^{(t + 1)}} - \nabla E \arrowvert_{w^{(t)}})}
 *    {\nabla E \arrowvert^T_{w^{(t + 1)}} \nabla E 
 *    \arrowvert_{w^{(t)}}}d^{(t)} 
 *  \f$
 *
 *  Because of the fact, that usually the target function is not quadratic
 *  and line search is not exact, unwished search directions can be
 *  calculated. So the algorithm can be restarted during optimization
 *  with the search direction set to the current negative gradient. <br>
 *
 *  For a more detailed description of this algorithm please refer to
 *  "Neural Networks for Pattern Recognition" by Christopher M. Bishop.
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
class CG : virtual public ModelInterface {

 public:

  //! Used to specify the type of line search algorithm
  //! that is used by the CG algorithm: <br>
  //! "Linmin" is the original line search algorithm,
  //! "DLinmin" is a line search algorithm that additionally uses
  //! derivative information and "Cblnsrch" is a cubic line 
  //! search.
  enum LineSearch { Dlinmin, Linmin, Cblnsrch };


  //! Initializes some internal variables used by the CG algorithm.
  void initCG( const Array<double> &in, const Array<double> &target, 
               LineSearch ls, int reset = -1, double ax = 0., double bx = 1., 
               double lambda = .25, bool verbose = false);


  //! Updates the weights of the network by using 
  //! the CG optimization algorithm and returns the minimum value by 
  //! reference.
  void cg(const Array<double> &in, const Array<double> &target, double& fret);


//===========================================================================
/*!
 *  \brief Updates the weights of the network by using 
 *         the CG optimization algorithm.
 *
 *  The given input patterns \em in and their corresponding target
 *  values \em target are used to calculate the error of the currently
 *  used model. The resulting error function is then minimized by the
 *  BFGS algorithm, the weights of the network are updated. <br>
 *  Notice: Method #initCG must be called before using the CG
 *  algorithm itself.
 *
 *      \param in         Input patterns for the currently used model.
 *      \param target     Target values corresponding to the input patterns.
 *      \return none
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
  void cg(const Array<double> &in, const Array<double> &target) 
  {
      double fret;
      cg(in, target, fret);
  };

  void reset();

 private:

   void dlinmin  (Array< double >& p,

                  const Array< double >& xi,

                  double&                fret,

                  const Array< double >& in,

                  const Array< double >& target);

   void linmin   (Array< double >& p,

                  const Array< double >& xi,

                  double&                fret,

                  const Array< double >& in,

                  const Array< double >& target);

   void cblnsrch (Array< double >& xold,

                  double fold,           // evaluation at current solution

                  Array< double >& g,    // current gradient

                  Array< double >& p,    // search direction

                  Array< double >& x,    // returns the new solution

                  double& fret,

                  const Array< double >&in,

                  const Array< double >&target); 





  



  LineSearch    lineSearchType;

  Array<double> CG_g, CG_h, CG_xi;

  unsigned      CG_n;



  double        CG_gg, CG_gam, CG_dgg;

  double        CG_fret;



  double        LS_lambda;

  double        LS_ax;

  double        LS_bx;



  unsigned      CG_count;

  unsigned      CG_reset;

  unsigned      CG_verbose;

};





#endif 


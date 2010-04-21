//===========================================================================
/*!
 *  \file CE.h
 *
 *  \brief Offers an (E)rror measure for (C)lassification tasks, that
 *         can be used for monitoring, but not for training.
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
 *      $RCSfile: CE.h,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:22 $
 *
 *  \par Changes:
 *      $Log: CE.h,v $
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.5  2002/05/16 13:24:23  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.4  2002/02/06 14:30:35  rudi
 *      Doxygen comments added.
 *
 *      Revision 1.3  2001/11/30 14:25:50  rudi
 *      Some doxygen comments added.
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

#ifndef CE_H
#define CE_H

#include <cmath>
#include "ReClaM/ModelInterface.h"

//===========================================================================
/*!
 *  \brief Offers an (E)rror measure for (C)lassification tasks, that
 *         can be used for monitoring, but not for training.
 * 
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
class CE : virtual public ModelInterface {
 public:

//===========================================================================
/*!
 *  \brief A classification error measure, similar to the "winner takes
 *         it all" scheme (ErrorMeasures::wta).
 *
 *  As "the winner takes it all" error measure (ErrorMeasures::wta)
 *  the output vectors \em model(in) are transformed into bitvectors and 
 *  compared to the target-bitvector \em out. In contrast to "wta", not 
 *  the output neuron
 *  with the best value is set to "1" and all other output neurons
 *  to zero, but all output neurons with a value greater than (max + min) / 2
 *  are set to max and all other output neurons 
 *  are set to min. 
 *
 *  \param in  Input vector for the model.
 *  \param out The target vector, that must be a bitvector.
 *  \param min Minimum value
 *  \param max Maximum value
 *  \return The classification error.
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
  double ce(const Array<double> &in, const Array<double> &out,
	    double min = 0, double max = 1.) {
    double ce = 0;
    if(in.ndim() == 1) {
      Array<double> output(out.dim(0));
      model(in, output);
      takeAll(output, min, max);
      if(out != output) ce++;
    } else {
      Array<double> output(out.dim(1));
      for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	model(in[pattern], output);
	takeAll(output, min, max);

	for(unsigned i = 0; i < output.nelem(); i++) {
	  if(out[pattern](i) != output(i)) {
	    ce++;
	    break;
	  }
	}
      }
      ce /= in.dim(0);
    }
    return ce;
  }


 private:

  // Transforms a linear vector into a bitvector by using
  // the threshold value of 0.5. All values greater than
  // 0.5 are set to 1 all other values to zero.
  //
  void takeAll(Array<double> &a, double min, double max) {
    for(unsigned i = 0; i < a.nelem(); i++) {
      if(a(i) > (max + min) / 2.) a(i) = max;
      else a(i) = min;
    }
  }
};

#endif

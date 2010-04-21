//===========================================================================
/*!
 *  \file BimodalBrownianProcess.h
 *
 *
 *  \author  Martin Kreutz
 *  \date    25.02.1999
 *
 *  \par Copyright (c) 1999, 2003:
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
 *      TestData
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: BimodalBrownianProcess.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: BimodalBrownianProcess.h,v $
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:07  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of Shark. This library is free software;
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


/*
 * This stochastic dynamical process decribes the motion of a particle that
 * is allowed to move freely in a double well potential:
 *
 *     V(x) = x^4 + 2 * x^2 + 1
 *
 * Due to the form of the potential this results in a bounded Brownian process.
 * Without external shocks the position of the particle would converge to one
 * of the two minima of the potential (-1 or 1).
 * The dynamic of the position is described by:
 *
 *     d^2 x     V'(x)           d x
 *     ----- = - ----- - alpha * --- + sigma * xi
 *     d t^2     mass            d t
 *
 * where 'mass' denotes the mass of the particle, 'alpha' the friction
 * coefficient and 'sigma' the standard deviation of the noise.
 * 'xi' is independent identically normally distributed with zero mean
 * and variance 1.
 *
 */


#ifndef __cplusplus
#error Must use C++.
#endif

#ifndef __BIMODALBROWNIANPROCESS_H
#define __BIMODALBROWNIANPROCESS_H

#ifdef __GNUC__
#pragma interface
#endif

#include "Rng/Normal.h"
#include "TestData/TimeSeries/RK4.h"

//===========================================================================

class BimodalBrownianProcess : public RK4
{
  public:
    BimodalBrownianProcess( double masspar  = 1,
			    double alphapar = 1,
			    double sigmapar = 1.5,
			    double steppar  = 0.1 );

    Array< double > operator ( ) ( );

  private:
    double mass, alpha, time, noise;
    Normal gauss;

    void derivative( const Array< double >&, Array< double >& ) const;

    Generator< Array< double > >* clone( ) const;
};

//===========================================================================

#endif /* !__BIMODALBROWNIANPROCESS_H */

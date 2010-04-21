//===========================================================================
/*!
 *  \file NoisyMackeyGlass.h
 *
 *
 *  \author  Martin Kreutz
 *  \date    16.09.1998
 *
 *  \par Copyright (c) 1998-2003:
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
 *      $RCSfile: NoisyMackeyGlass.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: NoisyMackeyGlass.h,v $
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:07  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of TestData. This library is free software;
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

/* NoisyMackeyGlass.h */
/* ======================================================================
 *
 *  File    :  NoisyMackeyGlass.h
 *  Created :  16.09.1998
 *
 *  Copyright (c) 1995,1998 Martin Kreutz
 *  EMail - Martin.Kreutz@neuroinformatik.ruhr-uni-bochum.de
 *
 *  Institut fuer Neuroinformatik
 *  Ruhr-Universitaet Bochum
 *  44780 Bochum
 *  
 *  Last update: 16.10.1998
 *
 * ----------------------------------------------------------------------
 *
 *  This file is part of the EALib. This library is free software;
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
 * ======================================================================
 */


#ifndef __cplusplus
#error Must use C++.
#endif

#ifndef __NOISYMACKEYGLASS_H
#define __NOISYMACKEYGLASS_H

#ifdef __GNUC__
#pragma interface
#endif

#include "Array/Array.h"
#include "Rng/Normal.h"
#include "TestData/TimeSeries/Generator.h"

//===========================================================================

class NoisyMackeyGlass : public Generator< double >
{
  public:
    NoisyMackeyGlass( double a, double b, double c, unsigned tau,
		      double sigma = 0, bool dyn = false);

    void   reset    ( );
    double operator ( ) ( );

    void   dynamicNoise ( ) { dynamic = true;  }
    void   additiveNoise( ) { dynamic = false; }

    void   setStdDev( double sigma ) { variance = sigma * sigma; }

  private:
    Normal   gauss;
    double   A, B, C;
    unsigned T, t;
    double   variance;
    bool     dynamic;
    Array< double > x;

    Generator< double >* clone( ) const;
};

//===========================================================================

#endif /* !__NOISYMACKEYGLASS_H */

//===========================================================================
/*!
 *  \file DiscreteMackeyGlass.h
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
 *      $RCSfile: DiscreteMackeyGlass.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: DiscreteMackeyGlass.h,v $
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


#ifndef __cplusplus
#error Must use C++.
#endif

#ifndef __DISCRETEMACKEYGLASS_H
#define __DISCRETEMACKEYGLASS_H

#ifdef __GNUC__
#pragma interface
#endif

#include "Array/Array.h"
#include "TestData/TimeSeries/Generator.h"

//===========================================================================

class DiscreteMackeyGlass : public Generator< double >
{
  public:
    DiscreteMackeyGlass( double   a   = 0.2,
			 double   b   = 10,
			 double   c   = 0.1,
			 unsigned tau = 30 );

    void   reset    ( );
    double operator ( ) ( );

  private:
    double   A, B, C; 
    unsigned T, t;
    Array< double > x;

    Generator< double >* clone( ) const;
};

//===========================================================================

#endif /* !__DISCRETEMACKEYGLASS_H */

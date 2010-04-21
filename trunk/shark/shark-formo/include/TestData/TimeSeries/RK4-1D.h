//===========================================================================
/*!
 *  \file RK4-1D.h
 *
 *
 *  \author  Martin Kreutz
 *  \date    15.10.1998
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
 *      $RCSfile: RK4-1D.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: RK4-1D.h,v $
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

#ifndef __RK4_1D_H
#define __RK4_1D_H

#ifdef __GNUC__
#pragma interface
#endif

#include "Array/Array.h"
#include "TestData/TimeSeries/Generator.h"

//===========================================================================

class RK4_1D : public Generator< double >
{
  public:
    RK4_1D( double stepsize );

    void reset( );
    double operator ( ) ( );

  protected:
    double   h;
    unsigned t;
    double y, dy, dym, dyt, yt;

    virtual double derivative( double ) const = 0;
};

//===========================================================================

#endif /* !__RK4_1D_H */

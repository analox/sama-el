//===========================================================================
/*!
 *  \file Ring.h
 *
 *
 *  \author  Martin Kreutz
 *  \date    20.08.1995
 *
 *  \par Copyright (c) 1995, 2003:
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
 *      $RCSfile: Ring.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: Ring.h,v $
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


#ifndef __cplusplus
#error Must use C++.
#endif

#ifndef __RING_H
#define __RING_H

#ifdef __GNUC__
#pragma interface
#endif

#include "Array/Array.h"
#include "Rng/Normal.h"
#include "Rng/Uniform.h"
#include "Mixture/RandomVector.h"

class Ring : public RandomVector< double >
{
  public:
    Ring( double xvar = 1. / 16,
	  double yvar = 1. / 16,
	  double xrad = 1,
	  double yrad = 1 );

    Ring( double xvar,
	  double yvar,
	  double xrad,
	  double yrad,
	  RNG& rng );

    double xvar    ( ) const         { return xStdDev * xStdDev; }
    double yvar    ( ) const         { return yStdDev * yStdDev; }
    double xrad    ( ) const         { return xRad;              }
    double yrad    ( ) const         { return yRad;              }
    void   xvar    ( double newVar ) { xStdDev = sqrt( newVar ); }
    void   yvar    ( double newVar ) { yStdDev = sqrt( newVar ); }
    void   xrad    ( double newRad ) { xRad    = newRad;         }
    void   yrad    ( double newRad ) { yRad    = newRad;         }

    Array< double > operator( )( );
    double          p( const Array< double >& ) const;

  private:
    Uniform uni;
    Normal  gauss;

    double  xStdDev, yStdDev;
    double  xRad,    yRad;

    void readFrom( std::istream& is );
    void writeTo ( std::ostream& os ) const;
};

#endif /* !__RING_H */

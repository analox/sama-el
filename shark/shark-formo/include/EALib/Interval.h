//===========================================================================
/*!
 *  \file Interval.h
 *
 *
 *  \author  Martin Kreutz
 *  \date    01.01.1995
 *
 *  \par Copyright (c) 1995-2003:
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
 *      EALib
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: Interval.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: Interval.h,v $
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of EALib. This library is free software;
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

#ifndef __INTERVAL_H
#define __INTERVAL_H

class Interval
{
  private:
    double   lower;
    double   upper;

  public:
    Interval ( )                      { lower = upper = 0;      }
    Interval ( double lo )            { lower = upper = lo;     }
    Interval ( double lo, double hi ) { lower = lo; upper = hi; }

    double   lowerBound ( ) const     { return lower;           }
    double   upperBound ( ) const     { return upper;           }
    double   width      ( ) const     { return upper - lower;   }
};

#endif /* !__INTERVAL_H */

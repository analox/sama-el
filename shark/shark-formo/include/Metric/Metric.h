//===========================================================================
/*!
 *  \file Metric.h
 *
 *
 *  \author  Axel W. Dietrich
 *  \date    2000-09-05
 *
 *  \par Copyright (c) 2000-2003:
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
 *      Metric
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: Metric.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: Metric.h,v $
 *      Revision 2.1  2004/03/04 12:55:28  shark-admin
 *      include file iostream.h exhchanged by iostream
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of Metric. This library is free software;
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

#ifndef __METRIC_H
#define __METRIC_H

#include <math.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

//===========================================================================

class Metric {
private:
  static inline int min( const int a, const int b, const int c );
  static inline int arrayLength( const int *array );

public:
  static int Levenshtein( const vector<bool> &a, const vector<bool> &b );
  static int Levenshtein( const string &a, const string &b );
  static int Levenshtein( const char *a, const char *b );
  static int Levenshtein( const int *a, const int *b );

  static int Hamming( const vector<bool> &a, const vector<bool> &b );
  static int Hamming( const string &a, const string &b );
  static int Hamming( const char *a, const char *b );
  static int Hamming( const int *a, const int *b );
};

#endif /* !__METRIC_H */

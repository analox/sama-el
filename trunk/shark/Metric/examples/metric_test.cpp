//===========================================================================
/*!
 *  \file metric_test.cpp
 *
 *
 *  \par Copyright (c) 2000,2003:
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
 *      $RCSfile: metric_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: metric_test.cpp,v $
 *      Revision 2.1  2004/03/04 13:34:11  shark-admin
 *      return value of 'int main()' added
 *
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
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


#include <iostream>
#include "Metric/Metric.h"

int main( )
{
  char *s1 = "0100010100";
  char *s2 = "1010111010";

  std::cout << std::endl << "Levenshtein distance" << std::endl;
  std::cout << "+ between " << s1 << " and " << s2 << " is " << Metric::Levenshtein( s1, s2 ) << std::endl;

  std::cout << std::endl;

  std::cout << "Hamming distance" << std::endl;
  std::cout << "+ between " << s1 << " and " << s2 << " is " << Metric::Hamming( s1, s2 ) << std::endl;

  std::cout << std::endl;

  return 1;
}

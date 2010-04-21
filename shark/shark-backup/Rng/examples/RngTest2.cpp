//===========================================================================
/*!
 *  \file RngTest2.cpp
 *
 *
 *  \par Copyright (c) 1999-2003:
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
 *      Rng
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: RngTest2.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: RngTest2.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:19  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:02  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of Rng. This library is free software;
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
#include "Rng/GlobalRng.h"
#include "Array/ArrayIo.h"

using namespace std;

int main(int argc, char **argv) {
  
	if (argc > 1) {
		Rng::seed(atoi(argv[1]));	
	}

	Array<unsigned> bin(10, 10);
	bin = 0;
  
	for(unsigned i = 0; i < 1000000; i++) 
	{
		unsigned a = Rng::discrete(0, 9);
		unsigned b = Rng::discrete(0, 9);
		bin(a, b) += 1; 
	}
	
	writeArray(bin, cout); 
	cout << endl;

	return 0;
}


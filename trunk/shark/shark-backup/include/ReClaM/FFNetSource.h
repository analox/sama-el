//===========================================================================
/*!
 *  \file FFNetSource.h
 *
 *  \author  M. Toussaint
 *  \date    2004
 *
 *  \par Copyright (c) 1999-2001:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-27974<BR>
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
 *      $RCSfile: FFNetSource.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/09 16:29:14 $
 *
 *  \par Changes:
 *      $Log: FFNetSource.h,v $
 *      Revision 2.1  2004/06/09 16:29:14  saviapbe
 *
 *      The construction "using namespace std;" was early. Corrected.
 *
 *      Revision 2.0  2004/05/27 15:12:24  saviapbe
 *
 *      Revision tad was reset to 2.x
 *
 *      Revision 1.1  2004/05/27 15:11:52  saviapbe
 *
 *      FFNetSource was added.
 *
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




#ifndef FFNET_SOURCE_H
#define FFNET_SOURCE_H

#include <string>
#include <stdio.h>
using namespace std;

extern "C" int sprintf(char *str, const char *format, ...);

void FFNetReplace(const string &str, string &newString, const string &token, const string &newToken1, int newToken2, const string &newToken3) {
  char buffer[256];

  sprintf(buffer, "%d", newToken2);
  string newToken = newToken1 + buffer + newToken3;

  string::size_type pos;
  newString = str;
  
  pos = newString.find (token, 0);
  while(pos != string::npos) {
    newString.replace(pos, token.length(), newToken);
    pos = newString.find (token, 0);
  }
}

void FFNetSource(ostream &os, unsigned inputDimension, unsigned outputDimension,
	   Array<int> &connectionMatrix, Array<double> &weightMatrix,
	   char *g, char *gOut, unsigned p = 4) {
  unsigned i, j;
  bool first;
  string str;

  unsigned numberOfNeurons = weightMatrix.dim(0);
  unsigned firstOutputNeuron = numberOfNeurons - outputDimension;
  unsigned bias = numberOfNeurons;

  os.precision(p);

  os << "void model(double *in, double *out) {\n";
  for(j = 0; j < inputDimension; j++) 
    os << "  double z" << j << " = in[" << j << "];\n";
  
  for(i = inputDimension; i < firstOutputNeuron; i++) {
    first = true;
    for(j = 0; j < i; j++) 
      if(connectionMatrix(i, j)) {
	if(first) {
	  os << "  double z" << i <<  " = z" << j << " * " << weightMatrix(i, j);
	  first = false;
	} else {
	  os << " + z" << j << " * " << weightMatrix(i, j);
	}
      }
    
    if(connectionMatrix(i, bias)) 
      if(first) {
	os  << "double z" << i << " = " << weightMatrix(i, bias); // bias
	first = false;
      } else {
	os << " + " << weightMatrix(i, bias);
      }

    if(!first) {
      os << ";\n";
      FFNetReplace(g,str,"#","z",i,"");
      os << "  z" << i << " = " << str << ";\n";
    }
  }
  for(i = firstOutputNeuron; i < numberOfNeurons; i++) {
    first = true;

    for(j = 0; j < i; j++)
      if(connectionMatrix(i, j)) {
	if(first) {
	  os << "  out[" << i-firstOutputNeuron << "] = " << "z" << j
	     << " * " << weightMatrix(i, j);
	  first = false;
	} else {
	  os << " + z" << j << " * " << weightMatrix(i, j);
	}
      }
    if(connectionMatrix(i, bias)) 
      if(first) {
	os << "  out[" << i-firstOutputNeuron << "] = " 
	   << weightMatrix(i, bias); // bias
	first = false;
      } else {
	os << " + " << weightMatrix(i, bias); // bias
      }
    if(!first) {
      os << ";\n";
      FFNetReplace(gOut,str,"#","out[",i-firstOutputNeuron,"]");
      os <<  "  out[" << i-firstOutputNeuron << "] = " << str << ";\n";
    }
  }
  os << "};" << endl;
}

#endif

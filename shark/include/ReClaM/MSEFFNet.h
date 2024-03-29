//===========================================================================
/*!
 *  \file MSEFFNet.h
 *
 *  \brief Offers the functions to create and to work with a
 *         feed-forward network combined with the mean squared error
 *         measure. This combination is created due to computational
 *         efficiency. 
 *         
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Copyright (c) 1999-2001:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
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
 *      $RCSfile: MSEFFNet.h,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:22 $
 *
 *  \par Changes:
 *      $Log: MSEFFNet.h,v $
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:26:03  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2002/02/06 14:46:23  rudi
 *      Doxygen comments added, memory leakage removed.
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

#ifndef MSEFFNET_H
#define MSEFFNET_H

#include "ReClaM/FFNet.h"

//===========================================================================
/*!
 *  \brief Offers the functions to create and to work with a
 *         feed-forward network combined with the mean squared error
 *         measure. This combination is created due to computational
 *         efficiency. 
 *         
 *  As with the "normal" feed forward networks, connections
 *  from neurons to the bias must be explicitly defined.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class MSEFFNet : public FFNet {

public:
  //! Creates an empty MSE feed-forward network with "in" 
  //! input neurons and "out" output neurons.
  MSEFFNet(const unsigned in = 0, const unsigned out = 0);


  //! Creates a MSE feed-forward network with "in" input neurons and
  //! "out" output neurons. Additionally, the array "cmat" determines
  //! the topology (i.e., number of neurons and their connections).
  MSEFFNet(const unsigned in, const unsigned out,
	    const Array<int>& cmat);


  //! Creates a MSE feed-forward network with "in" input neurons and
  //! "out" output neurons. Additionally, the arrays "cmat" and
  //! "wmat" determine the topology (i.e., number of neurons and their
  //! connections) as well as the connection weights.
  MSEFFNet(const unsigned in, const unsigned out,
	    const Array<int>& cmat, const Array<double>& wmat);
  
  //! Creates a MSE feed-forward network by reading the necessary 
  //! information from a file named "filename".
  MSEFFNet(const std::string &);


  //! Destructs an MSE Feed Forward Network object.
  virtual ~MSEFFNet() {};


  //! Calculates the mean squared error of the model output, compared to 
  //! the target values. 
  double error(const Array<double> &input, const Array<double> &target);


  //! Calculates the mean squared Error of the model output, 
  //! compared to the target values, and the derivative of the 
  //! mean squared error with respect to the weights.
  double derror(const Array<double> &input, const Array<double> &target,
		bool returnError = true);


  //! Method is not working yet. Will always exit with failure.
  void dmodel(const Array<double> &input);


  //! Method is not working yet. Will always exit with failure.
  void dmodel(const Array<double> &input, Array<double> &output);

 protected:
  void resize() { FFNet::resize(); d.resize(numberOfNeurons); }
};


//! This is done for backwards compatibility. Now, that
//! bias connections must be explicitly defined, there
//! is no (declaration) difference between a "MSEFFNet"
//! and a "MSEBFFNet".
typedef  MSEFFNet MSEBFFNet;

#endif









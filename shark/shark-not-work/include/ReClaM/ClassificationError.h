//===========================================================================
/*!
*  \file ClassificationError.h
*
*  \brief Compute the fraction of classification errors
*
*  \author  T. Glasmachers
*  \date    2006
*
*  \par Copyright (c) 1999-2006:
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
 *      $RCSfile: ClassificationError.h,v $<BR>
*
*  \par Changes:
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

#ifndef _ClassificationError_H_
#define _ClassificationError_H_


#include <ReClaM/ErrorFunction.h>


//! The ClassificationError class returns the number of
//! classification errors. This measure is not differentiable.
class ClassificationError : public ErrorFunction
{
public:
	//! Constructor
	ClassificationError()
	{
	}

	//! Destructor
	~ClassificationError()
	{
	}


	//! Computation of the fraction of wrongly classified examples.
	double error(Model& model, const Array<double>& input, const Array<double>& target)
	{
		int wrong = 0;
		int i, ic = input.dim(0);
		Array<double> output;
		for (i=0; i<ic; i++)
		{
			model.model(input[i], output);
			if (output(0) * target(i) <= 0.0) wrong++;
		}
		return ((double)wrong) / ic;
	}
};


//! The ClassificationError class returns the number of
//! classification errors, rescaled by the class magnitudes.
//! For unbalanced datasets, it is in most cases preferable
//! compared to the #ClassificationError.
class BalancedClassificationError : public ErrorFunction
{
public:
	//! Constructor
	BalancedClassificationError()
	{
	}

	//! Destructor
	~BalancedClassificationError()
	{
	}


	//! Computation of the fraction of wrongly classified examples,
	//! rescaled as follows:
	//! let p be the number of positive class examples and let
	//! n be the number of negative class examples. Let a be the
	//! number of positive examples classified negative (type 1 error)
	//! and let b be the number of negative examples classifies positive
	//! (type 2 error). Then, the balanced classification error is the
	//! error rate
	//! \f[
	//!       \frac{1}{2} \left( \frac{a}{p} + \frac{b}{n} \right) .
	//! \f]
	//! The expectation value of the balanced classification error is
	//! invariant under the class magnitude, that is the fraction of
	//! positive and negative class examples. For example, if all
	//! positive examples are doubled in the test set, the error rate
	//! remains unchanged.
	double error(Model& model, const Array<double>& input, const Array<double>& target)
	{
		SIZE_CHECK(input.ndim() == 2);
		SIZE_CHECK(target.ndim() == 1);

		int p = 0;
		int n = 0;
		int a = 0;
		int b = 0;
		int i, ic = input.dim(0);
		SIZE_CHECK(target.dim(0) == ic);

		for (i=0; i<ic; i++)
		{
			Array<double> output;
			model.model(input[i], output);
			if (target(i) > 0.0)
			{
				p++;
				if (output(0) <= 0.0) a++;
			}
			else
			{
				n++;
				if (output(0) > 0.0) b++;
			}
		}

		return 0.5 * ((((double)a) / p) + (((double)b) / n));
	}
};


#endif

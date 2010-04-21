//===========================================================================
/*!
*  \file optGrid.h
*
*  \brief optimization by grid or point set search
*
*  This file provides a collection of quite simple optimizers.
*  It provides a basic grid search, a nested grid search and
*  a search on a predefined set of points.
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
*      $RCSfile: optGrid.h,v $<BR>
*
*  \par Changes:
*      $Log: optGrid.h,v $
*      Revision 2.7  2006/04/20 16:25:25  glasmtbl
*      *** empty log message ***
*
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

#ifndef _optGrid_H_
#define _optGrid_H_


#include <ReClaM/Optimizer.h>
#include <iostream>


//!
//! \brief Optimize by trying out a grid of configurations
//!
//! \par
//! The optGrid class allows for the definition of a grid in
//! parameter space. It does a simple one-step optimization over
//! the grid by trying out every possible parameter combination.
//! Please note that the computation effort grows exponentially
//! with the number of parameters.
//!
//! \par
//! If you only want to try a subset of the grid, consider using
//! the #optPointSearch class instead.
//! A more sophisticated (less exhaustive) grid search variant is
//! available with the #optGridNested class.
//!
class optGrid
{
public:
	//! Constructor
	optGrid()
	{
	}

	//! Destructor
	~optGrid()
	{
	}


	//! basic initialization with default parameters
	void init(Model& model)
	{
		throw "[optGrid::init] A default initialization is impossible for the grid optimizer.";
	}

	//! uniform initialization for all parameters
	//! \param  params  number of model parameters to optimize
	//! \param  min     smallest parameter value
	//! \param  max     largest parameter value
	//! \param  nodes   total number of values in the interval
	void init(int params, double min, double max, int nodes)
	{
		numberOfValues.resize(params, false);
		nodeValues.resize(params, nodes, false);
		int i, j;
		for (j=0; j<params; j++)
		{
			for (i=0; i<nodes; i++)
			{
				nodeValues(j, i) = min + i * (max - min) / (nodes - 1.0);
			}
		}
	}

	//! individual definition for every parameter
	//! \param  params  number of model parameters to optimize
	//! \param  min     smallest value for every parameter
	//! \param  max     largest value for every parameter
	//! \param  nodes   total number of values for every parameter
	void init(int params, const Array<double>& min, const Array<double>& max, const Array<int>& nodes)
	{
		int i, j;
		int nmax = 0;
		for (i=0; i<params; i++) if (nodes(i) > nmax) nmax = nodes(i);
		numberOfValues.resize(params, false);
		nodeValues.resize(params, nmax, false);
		for (j=0; j<params; j++)
		{
			for (i=0; i<nodes(j); i++)
			{
				nodeValues(j, i) = min(j) + i * (max(j) - min(j)) / (nodes(j) - 1.0);
			}
		}
	}

	//! uniform definition of the values to test for all parameters
	//! \param  params  number of model parameters to optimize
	//! \param  values  values used for every coordinate
	void init(int params, const Array<double>& values)
	{
		int d = values.dim(0);
		int i, j;
		numberOfValues.resize(params, false);
		numberOfValues = d;
		nodeValues.resize(params, d, false);
		for (j=0; j<params; j++)
		{
			for (i=0; i<d; i++)
			{
				nodeValues(j, i) = values(i);
			}
		}
	}

	//! individual definition for every parameter
	//! \param  params  number of model parameters to optimize
	//! \param  nodes   total number of values for every parameter
	//! \param  values  values used. The first dimension is the parameter, the second dimension is the node. Unused entries are ignored.
	void init(int params, const Array<int>& nodes, const Array<double>& values)
	{
		if (params != nodes.dim(0) || params != values.dim(0))
			throw "[ConcatenatedModel::AppendModel] Dimension conflict";

		numberOfValues = nodes;
		nodeValues = values;
	}

	//! Please note that for the grid search optimizer it does
	//! not make sense to call #optimize more than once, as the
	//! solution does not improve iteratively.
	double optimize(Model& model, ErrorFunction& errorfunction, const Array<double>& input, const Array<double>& target)
	{
		int params = numberOfValues.dim(0);
		Array<int> index(params);
		index = 0;
		int i;
		double e;
		double best = 1e100;
		Array<double> best_param(params);

		// loop through all grid points
		while (true)
		{
			// define the parameters
			for (i=0; i<params; i++) model.setParameter(i, nodeValues(i, index(i)));

			// evaluate the model
			if (model.isFeasible())
			{
				e = errorfunction.error(model, input, target);
				if (e < best)
				{
					best = e;
					for (i=0; i<params; i++) best_param(i) = model.getParameter(i);
				}
			}

			// next index
			for (i=params-1; i>=0; i--)
			{
				index(i)++;
				if (index(i) == numberOfValues(i))
				{
					index(i) = 0;
				}
				else break;
			}
			if (i == -1) break;
		}

		// write the best parameters into the model
		for (i=0; i<params; i++) model.setParameter(i, best_param(i));
	}

protected:
	//! The array holds the number of grid values for every parameter axis.
	Array<int> numberOfValues;

	//! The array columns contain the grid values for the corresponding parameter axis.
	//! As all columns have the same size, some values may be meaningless.
	Array<double> nodeValues;
};


//!
//! \brief Nested grid search
//!
//! \par
//! The optGridNested class is an iterative optimizer,
//! doing one grid search in every iteration. In every
//! iteration, it halves the grid extent doubling the
//! resolution in every coordinate.
//!
//! \par
//! Although nested grid search is much less exhaustive
//! than standard grid search, it still suffers from
//! exponential time and memory complexity in the number
//! of variables optimized. Therefore, if the number of
//! variables is larger than 2 or 3, consider using the
//! #optCMA class instead.
//!
class optGridNested : public Optimizer
{
public:
	//! Constructor
	optGridNested()
	{
	}

	//! Destructor
	~optGridNested()
	{
	}


	//! There is no useful default initialization for this optimizer.
	//! Thus, this member throws an exception.
	void init(Model& model)
	{
		throw "[optGridNested::init] A default initialization is impossible for the grid optimizer.";
	}

	//!
	//! \brief Initialization of the nested grid search.
	//!
	//! \par
	//! The min and max arrays define ranges for every parameter to optimize.
	//! These ranges are strict, that is, the algorithm will not try values
	//! beyond the range, even if is finds a boundary minimum.
	//!
	//! \param  model  #Model to optimize
	//! \param  min    lower end of the parameter range
	//! \param  max    upper end of the parameter range
	void init(Model& model, const Array<double>& min, const Array<double>& max)
	{
		SIZE_CHECK(min.ndim() == 1);
		SIZE_CHECK(max.ndim() == 1);

		int d, dc = min.dim(0);
		SIZE_CHECK(max.dim(0) == dc);
		RANGE_CHECK(model.getParameterDimension() >= dc);

		minimum = min;
		maximum = max;

		step.resize(dc, false);
		std::vector<unsigned int> D(dc);
		for (d=0; d<dc; d++)
		{
			model.setParameter(d, 0.5 * (min(d) + max(d)));
			step(d) = 0.25 * (max(d) - min(d));
			D[d] = 5;
		}

		landscape.resize(D, false);
		landscape = 1e100;
	}

	//! Every call of the optimization member computes the
	//! error landscape on the current grid. It picks the
	//! best error value and zooms into the error landscape
	//! by a factor of 2 reusing already computed grid points.
	double optimize(Model& model, ErrorFunction& errorfunction, const Array<double>& input, const Array<double>& target)
	{
		SIZE_CHECK(step.ndim() == 1);
		SIZE_CHECK(minimum.ndim() == 1);
		SIZE_CHECK(maximum.ndim() == 1);

		int d, dc = step.dim(0);
		bool bNew;
		double value;
		double e, best = 1e99;
		std::vector<unsigned int> index(dc);
		std::vector<unsigned int> best_index(dc);
		Array<double> best_param(dc);
		Array<double> old_param(dc);

		SIZE_CHECK(minimum.dim(0) == dc);
		SIZE_CHECK(maximum.dim(0) == dc);
		RANGE_CHECK(model.getParameterDimension() >= dc);

		// initialize variables
		for (d=0; d<dc; d++)
		{
			value = model.getParameter(d);
			old_param(d) = value;
			best_param(d) = value;
			index[d] = 0;
			best_index[d] = 2;
		}

		// loop through the grid
		while (true)
		{
			// compute the grid point,
			// define it as the model parameters
			// and check whether the computation
			// has to be done at all
			bNew = true;
			e = landscape(index);

			if (e < 1e100) bNew = false;
			else
			{
				// set the parameters
				for (d=0; d<dc; d++)
				{
					value = old_param(d) + (index[d] - 2.0) * step(d);
					if (value < minimum(d) || value > maximum(d)) bNew = false;
					model.setParameter(d, value);
				}
			}

			// evaluate the grid point
			if (bNew)
			{
				if (model.isFeasible())
				{
					e = errorfunction.error(model, input, target);
					landscape(index) = e;
				}
			}

			// remember the best solution
			if (e < best)
			{
				best = e;
				for (d=0; d<dc; d++)
				{
					best_index[d] = index[d];
					best_param[d] = old_param(d) + (index[d] - 2.0) * step(d);
				}
			}

			// move to the next grid point
			for (d=0; d<dc; d++)
			{
				index[d]++;
				if (index[d] <= 4) break;
				index[d] = 0;
			}
			if (d == dc) break;
		}

		// zoom into the error landscape array
		Array<double> zoomed = landscape;
		std::vector<unsigned int> zoomed_index(dc);
		for (d=0; d<dc; d++) zoomed_index[d] = 0;
		// loop through the grid
		while (true)
		{
			// if appropriate, copy the error landscape point
			for (d=0; d<dc; d++)
			{
				if ((zoomed_index[d] & 1) == 1) break;
				index[d] = best_index[d] + (zoomed_index[d] / 2) - 1;
				if (index[d] > 4) break;
			}
			if (d == dc) zoomed(zoomed_index) = landscape(index);
			else zoomed(zoomed_index) = 1e100;

			// move to the next grid point
			for (d=0; d<dc; d++)
			{
				zoomed_index[d]++;
				if (zoomed_index[d] <= 4) break;
				zoomed_index[d] = 0;
			}
			if (d == dc) break;
		}
		landscape = zoomed;

		// decrease the step sizes
		for (d=0; d<dc; d++) step(d) *= 0.5;

		// load the best parameter configuration into the model
		for (d=0; d<dc; d++) model.setParameter(d, best_param(d));

		// return the lowest error value archieved
		return best;
	}

protected:
	//! minimum parameter value to check
	Array<double> minimum;

	//! maximum parameter value to check
	Array<double> maximum;

	//! current step size for every parameter
	Array<double> step;

	//! \brief error landscape
	//!
	//! \par
	//! The landscape array has as many dimensions as
	//! there are parameters to optimize. In every array
	//! dimension there are 5 entries. During every grid
	//! search iteration, the error is computed for all
	//! grid points not seen during previous iterations.
	//!
	//! \par
	//! Let N denote the number of parameters to optimize.
	//! To compute the error ladndscape at the current
	//! zoom level, the algorithm has to do \f$ 5^N \f$
	//! error function evaluation in the first iteration,
	//! and up to \f$ 5^N - 3^N \f$ evaluations in
	//! subsequent iterations.
	//!
	//! \par
	//! The grid is always centered around the best
	//! solution currently known. If this solution is
	//! located at the boundary, the landscape may exceed
	//! the parameter range defined #minimum and #maximum.
	//! These invalid landscape values are not used and
	//! are always set to 1e100, indicating non-optimality.
	Array<double> landscape;
};


//!
//! \brief Optimize by trying out predefined configurations
//!
//! \par
//! The optPointSearch class is similair to the #optGrid class
//! by the property that it optimizes a model in a single pass
//! just trying out a predefined number of parameter configurations.
//! The main difference is that every parameter configuration has
//! to be explicitly defined. It is not possible to define a set
//! of values for every axis; refer to #optGrid for this purpose.
//! Thus, the optPointSearch class allows for more flexibility.
//!
class optPointSearch : public Optimizer
{
public:
	//! Constructor
	optPointSearch()
	{
	}

	//! Destructor
	~optPointSearch()
	{
	}


	//! basic initialization with default parameters
	void init(Model& model)
	{
		throw "[optPointSearch::init] A default initialization is impossible for the point search optimizer.";
	}

	//! Initialization of the search points.
	//! \param  values  two-dimensional array; every column contains one parameter configuration.
	void init(Array<double>& values)
	{
		nodes = values;
	}

	//! Please note that for the point search optimizer it does
	//! not make sense to call #optimize more than once, as the
	//! solution does not improve iteratively.
	double optimize(Model& model, ErrorFunction& errorfunction, const Array<double>& input, const Array<double>& target)
	{
		int t, tc = nodes.dim(0);
		int p, pc = nodes.dim(1);
		double e;
		double best = 1e100;
		int best_index = -1;

		// loop through all points
		for (t=0; t<tc; t++)
		{
			// define the parameters
			for (p=0; p<pc; p++) model.setParameter(p, nodes(t, p));

			// evaluate the model
			if (model.isFeasible())
			{
				e = errorfunction.error(model, input, target);
				if (e < best)
				{
					best = e;
					best_index = t;
				}
			}
		}

		// write the best parameters into the model
		for (p=0; p<pc; p++) model.setParameter(p, nodes(best_index, p));
	}

protected:
	//! The array holds one parameter configuration in every column.
	Array<double> nodes;
};


#endif

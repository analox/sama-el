//===========================================================================
/*!
 *  \file CrossValidation.h
 *
 *  \brief Cross Validation
 *
 *
 *  The cross validation procedure is designed to adapt so-called
 *  hyperparameters of a model, which is based upon another model.
 *  In every hyperparameter evaluation step, the base model is
 *  trained. This inner training procedure depends on the
 *  hyperparameters. Thus, the hyperparameters can be evaluated by
 *  simply evaluating the trained base model. To avoid empirical
 *  risk minimization, the cross validation procedure splits the
 *  available data into training and validation subsets, such that
 *  all data points appear in training and validation subsets
 *  equally often.
 *
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

#ifndef _CrossValidation_H_
#define _CrossValidation_H_


#include "ReClaM/Model.h"
#include "ReClaM/ErrorFunction.h"
#include "ReClaM/Optimizer.h"
#include "Rng/GlobalRng.h"
#include <vector>


//! The Partitioning class defined a partitioning of
//! a set of training points and labels as it is required
//! for a cross validation procedure. That is, the data
//! are split into N subsets, usually of comparable size.
//! Then, the data are combined into training and validation
//! sets in N different ways, called partitions. For every
//! partition, one of the N subsets is used for validation,
//! while the union of all other subsets it used for
//! training.
class Partitioning
{
public:
	//!
	//! \brief Constructor
	//!
	//! The Constructor constructs the subsets.
	//! Every subset contains (approximately) the same
	//! number of elements. For every partition, all
	//! but one subset form the training set, while the
	//! remaining one is used for validation.
	//!
	//! \param numberOfPartitions  number of partitions to create
	//! \param input               input data to separate into subsets
	//! \param target              labels corresponding to the input data
	//!
	Partitioning(int numberOfPartitions, const Array<double>& input, const Array<double>& target)
	{
		partitions = numberOfPartitions;
		int ec = input.dim(0);
		int dimension = input.dim(1);
		int i, e, d;

		SIZE_CHECK(input.ndim() == 2);
		SIZE_CHECK(target.ndim() == 2);
		SIZE_CHECK(target.dim(0) == ec);

		part_index.resize(ec, false);
		part_train_input.resize(partitions, false);
		part_train_target.resize(partitions, false);
		part_validation_input.resize(partitions, false);
		part_validation_target.resize(partitions, false);

		int l_plus = 0;
		int l_minus = 0;
		for (e=0; e<ec; e++) if (target(e) > 0.0) l_plus++; else l_minus++;

		Array<int> t(partitions);
		Array<int> v_plus(partitions);
		Array<int> v_minus(partitions);
		int nn_plus = l_plus / partitions;
		int c_plus = l_plus - nn_plus * partitions;
		int nn_minus = l_minus / partitions;
		int c_minus = l_minus - nn_minus * partitions;
		for (i=0; i<partitions; i++)
		{
			v_plus(i) = nn_plus;
			v_minus(i) = nn_minus;
			if (i < c_plus) v_plus(i)++;
			if (i < c_minus) v_minus(i)++;
			t(i) = ec - v_plus(i) - v_minus(i);
			part_train_input[i] = new Array<double>(t(i), dimension);
			part_train_target[i] = new Array<double>(t(i));
			part_validation_input[i] = new Array<double>(v_plus(i) + v_minus(i), dimension);
			part_validation_target[i] = new Array<double>(v_plus(i) + v_minus(i));
		}
		int c;
		for (e=0; e<ec; e++)
		{
			if (target(e) > 0.0)
			{
				do
				{
					c = Rng::discrete(0, partitions - 1);
				}
				while (v_plus(c) == 0);
				v_plus(c)--;
			}
			else
			{
				do
				{
					c = Rng::discrete(0, partitions - 1);
				}
				while (v_minus(c) == 0);
				v_minus(c)--;
			}
			part_index(e) = c;

			for (d=0; d<dimension; d++) part_validation_input[c]->operator () (v_plus(c) + v_minus(c), d) = input(e, d);
			part_validation_target[c]->operator () (v_plus(c) + v_minus(c)) = target(e);

			for (i=0; i<partitions; i++)
			{
				if (i == c) continue;
				t(i)--;
				for (d=0; d<dimension; d++) part_train_input[i]->operator () (t(i), d) = input(e, d);
				part_train_target[i]->operator () (t(i)) = target(e);
			}
		}
	}

	//! Destructor
	~Partitioning()
	{
		int i;
		for (i=0; i<partitions; i++)
		{
			delete part_train_input[i];
			delete part_train_target[i];
			delete part_validation_input[i];
			delete part_validation_target[i];
		}
	}


	//! Return the number of partitions in the partitioning.
	inline int getNumberOfPartitions()
	{
		return partitions;
	}

	//! Return the index-th training point array
	inline const Array<double>& train_input(int index)
	{
		return *(part_train_input[index]);
	}

	//! Return the index-th training label array
	inline const Array<double>& train_target(int index)
	{
		return *(part_train_target[index]);
	}

	//! Return the index-th validation point array
	inline const Array<double>& validation_input(int index)
	{
		return *(part_validation_input[index]);
	}

	//! Return the index-th validation label array
	inline const Array<double>& validation_target(int index)
	{
		return *(part_validation_target[index]);
	}

	//! Return the validation set index for a given example
	inline int getIndex(int example)
	{
		return part_index(example);
	}

protected:
	//! number of partitions of the data
	int partitions;

	//! index array: which example belongs to which validation set?
	Array<int> part_index;

	//! training points
	std::vector< Array<double>* > part_train_input;

	//! training labels
	std::vector< Array<double>* > part_train_target;

	//! validation points
	std::vector< Array<double>* > part_validation_input;

	//! validation labels
	std::vector< Array<double>* > part_validation_target;
};


//!
//! \brief Collection of sub-models for cross validation
//!
//! \par
//! The CVModel class is based upon a set of models. The idea is
//! to have one independent model for every partition of the data
//! during the cross validation procedure. The class simply
//! collects these base models and synchronizes its parameters,
//! as long as they are accessed via the CVModel class.
//!
//! \par
//! In principle, it is not clear which base model to use for
//! prediction on unseen data. However, it makes sence to use the
//! base models on the whole cross validation data set, that is
//! as well on the subsets used for training, as well as on the
//! subsets for validation.
class CVModel : public Model
{
public:
	//! Constructor
	//!
	//! \param  models  Array of identical base models. There must be as many models in the array as there are cross validation partitions.
	CVModel(Array<Model*>& models)
	{
		SIZE_CHECK(models.ndim() == 1);
		SIZE_CHECK(models.dim(0) > 1);

		int p, pc = models(0)->getParameterDimension();

#ifdef DEBUG
		int i, ic = models.dim(0);
		for (i=1; i<ic; i++) SIZE_CHECK(models.getParameterDimension() == pc);
#endif

		baseModel = models;
		baseModelIndex = 0;

		parameter.resize(pc);
		for (p=0; p<pc; p++) parameter(p) = models(0)->getParameter(p);
	}

	//! Destructor
	~CVModel()
	{
	}


	//! Modifies a specific model parameter.
	void setParameter(unsigned int index, double value)
	{
		Model::setParameter(index, value);

		int i, ic = baseModel.dim(0);
		for (i=0; i<ic; i++) baseModel(i)->setParameter(index, value);
	}

	//! Set the currently used base model.
	void setBaseModel(int index)
	{
		RANGE_CHECK(index >= 0 && index < baseModel.dim(0));
		baseModelIndex = index;
	}

	//! Return a reference to the current base model.
	inline Model& getBaseModel()
	{
		return *baseModel(baseModelIndex);
	}

	//! The base model is used for model computations.
	void model(const Array<double>& input, Array<double>& output)
	{
		baseModel(baseModelIndex)->model(input, output);
	}

	//! The base model is used for model computations.
	void modelDerivative(const Array<double>& input, Array<double>& derivative)
	{
		baseModel(baseModelIndex)->modelDerivative(input, derivative);
	}

	//! The base model is used for model computations.
	void modelDerivative(const Array<double>& input, Array<double>& output, Array<double>& derivative)
	{
		baseModel(baseModelIndex)->modelDerivative(input, output, derivative);
	}

	//! Check whether the parameters define a valid model
	bool isFeasible()
	{
		return getBaseModel().isFeasible();
	}

protected:
	//! Model used for every data partition. Use the members
	//! #SetModel and #StoreModel to handle parameter exchange.
	Array<Model*> baseModel;

	//! Index of the baseModel to use
	int baseModelIndex;
};


//!
//! \brief #ErrorFunction based on a cross validation procedure
//!
//! \par
//! The CVError class computes the mean error over partitions,
//! that is, the cross validation error. It trains the sub-models
//! defined by the #CVModel object with a given #Optimizer on a
//! given #ErrorFunction using the training parts of the partitions.
//! For the mean error computation, it uses the varidation part of
//! the partitions, which contain every example exactly once.
class CVError : public ErrorFunction
{
public:
	//! Constructor
	//!
	//! \param  part       #Partitioning defining subsets
	//! \param  error      #ErrorFunction used for the subset tasks
	//! \param  optimizer  #Optimizer used for the subset tasks
	//! \param  iter       number of optimization iterations for the subset tasks
	CVError(Partitioning& part, ErrorFunction& error, Optimizer& optimizer, int iter)
	: partitioning(part)
	, baseError(error)
	, baseOptimizer(optimizer)
	{
		iterations = iter;
	}

	//! Destructor
	~CVError()
	{
	}


	//! Compute the cross validation error defined as the mean
	//! error over the subsets. The model parameters remain
	//! unchainged, as they are restored after every training
	//! procedure.
	//!
	//! \param  model  The model parameter has to be a reference to a #CVModel object.
	double error(Model& model, const Array<double>& input, const Array<double>& target)
	{
		CVModel* pCVM = dynamic_cast<CVModel*>(&model);
		if (pCVM == NULL) throw "[CVError::error] invalid model";

		SIZE_CHECK(input.ndim() == 2);
		SIZE_CHECK(target.ndim() == 1);
		SIZE_CHECK(target.dim(0) == input.dim(0));

		double err, ret = 0.0;
		int t, it, parts = partitioning.getNumberOfPartitions();
		int p, pc = model.getParameterDimension();
		Array<double> train_input;
		Array<double> train_target;
		Array<double> validation_input;
		Array<double> validation_target;
		Array<double> initial_param(pc);

		// loop through the partitions
		// std::cout << "[" << std::flush;
		for (t=0; t<parts; t++)
		{
			// std::cout << "<" << std::flush;
			// activate the t-th submodel
			pCVM->setBaseModel(t);

			// do the training
			baseOptimizer.init(pCVM->getBaseModel());
			for (it=0; it<iterations; it++)
			{
				baseOptimizer.optimize(pCVM->getBaseModel(), baseError, partitioning.train_input(t), partitioning.train_target(t));
			}

			// compute the validation error
			err = baseError.error(pCVM->getBaseModel(), partitioning.validation_input(t), partitioning.validation_target(t));
			ret += err;
			// std::cout << err << ">" << std::flush;
		}
		// std::cout << "]" << std::flush;

		// return the mean error
		return ret / parts;
	}

	double errorDerivative(Model& model, const Array<double>& input, const Array<double>& target, Array<double>& derivative)
	{
		CVModel* pCVM = dynamic_cast<CVModel*>(&model);
		if (pCVM == NULL) throw "[CVError::errorDerivative] invalid model";

		SIZE_CHECK(input.ndim() == 2);
		SIZE_CHECK(target.ndim() == 1);
		SIZE_CHECK(target.dim(0) == input.dim(0));

		double err, ret = 0.0;
		int t, it, parts = partitioning.getNumberOfPartitions();
		int p, pc = model.getParameterDimension();
		Array<double> train_input;
		Array<double> train_target;
		Array<double> validation_input;
		Array<double> validation_target;
		Array<double> initial_param(pc);
		Array<double> innerDerivative;

		derivative.resize(pc, false);
		derivative = 0.0;

		// loop through the partitions
		for (t=0; t<parts; t++)
		{
			// activate the t-th submodel
			pCVM->setBaseModel(t);

			// do the training
			baseOptimizer.init(pCVM->getBaseModel());
			for (it=0; it<iterations; it++)
			{
				baseOptimizer.optimize(pCVM->getBaseModel(), baseError, partitioning.train_input(t), partitioning.train_target(t));
			}

			// compute the validation error
			err = baseError.errorDerivative(pCVM->getBaseModel(), partitioning.validation_input(t), partitioning.validation_target(t), innerDerivative);
			for (p=0; p<pc; p++) derivative(p) += innerDerivative(p);
			ret += err;
		}

		// return the mean derivative
		for (p=0; p<pc; p++) derivative(p) = derivative(p) / parts;

		// return the mean error
		return ret / parts;
	}

protected:
	//! Partitioning upon which the cross validation procedure is based
	Partitioning& partitioning;

	//! ErrorFunction to use for the single cross validation tasks
	ErrorFunction& baseError;

	//! Optimizer to use for the single cross validation tasks
	Optimizer& baseOptimizer;

	//! Number of iterations to perform for the optimization of the single cross validation tasks
	int iterations;
};


#endif

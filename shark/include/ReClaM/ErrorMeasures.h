//===========================================================================
/*!
 *  \file ErrorMeasures.h
 *
 *  \brief Various error measures (not for training!).
 *         
 *  This file offers a bunch of error measures that can be used for 
 *  calculating errors in a different way for validation, but not for the
 *  training process.
 *
 *  \author  C. Igel
 *  \date    1999
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
 *      $RCSfile: ErrorMeasures.h,v $<BR>
 *      $Revision: 2.2 $<BR>
 *      $Date: 2005/08/03 12:02:33 $
 *
 *  \par Changes:
 *      $Log: ErrorMeasures.h,v $
 *      Revision 2.2  2005/08/03 12:02:33  christian_igel
 *      changed limits/values etc.
 *
 *      Revision 2.1  2004/06/01 15:41:45  saviapbe
 *
 *      The bugs in the doxygen's documentation were removed.
 *
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.5  2002/05/16 13:24:39  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.4  2001/11/30 14:45:13  rudi
 *      obsolete command removed
 *      doxygen comments added.
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

#ifndef ERROR_MEASURES_H
#define ERROR_MEASURES_H

#include <cmath>

#ifndef __SOLARIS__
    #include <limits>
#else
    #include <values.h>
#endif

#include "ReClaM/ModelInterface.h"
#include "Array/ArrayOp.h"

//===========================================================================
/*!
 *  \brief Various error measures (not for training!).
 *         
 *  This file offers a bunch of error measures that can be used for calculating
 *  errors in a different way for validation, but not for the
 *  training process.
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
class ErrorMeasures : virtual public ModelInterface {
 public:

//===========================================================================
/*!
 *  \brief Calculates the error percentage based on the mean squared error.
 *
 *  Measures the euklidian distance between the model output \em model(in),
 *  calculated from the input vector \em in, and the target vector \em out. 
 *  Besides normalizing the squared error to the number \em P of patterns,
 *  it is also divided with the number of outputs \em N and the range
 *  of the target vectors. The resulting value is used to determine
 *  the percentage of errors.
 *  \f[
 *      E = \frac{100}{NP(max\{out_{ip}\} - min\{out_{ip}\})^{2}}
 *      \sum_{p=1}^P\sum_{i=1}^N(model(in)_{ip} - 
 *      out_{ip})^{2}
 *  \f]
 *
 *      \param  in Input vector for the model.
 *      \param  out Target vector.
 *      \return The error percentage \em E.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
 double errorPercentage(const Array<double> &in, const Array<double> &out) 
{
    double se = 0;
#ifndef __SOLARIS__
    double outmax = std::numeric_limits< double >::min( );
    double outmin = std::numeric_limits< double >::max( );
#else
    double outmax = -MAXDOUBLE;
    double outmin = MAXDOUBLE;
#endif
    if(in.ndim() == 1) {
      Array<double> output(out.dim(0));
      model(in, output);
      for(unsigned c = 0; c < out.dim(0); c++) {
	if(out(c) > outmax) outmax = out(c);
	if(out(c) < outmin) outmin = out(c);
	se += (out(c) - output(c)) * (out(c) - output(c));
      }
    } else {
      Array<double> output(out.dim(1));
      for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	model(in[pattern], output);
	for(unsigned c = 0; c < out.dim(1); c++) {
	  if(out(pattern, c) > outmax) outmax = out(pattern, c);
	  if(out(pattern, c) < outmin) outmin = out(pattern, c);
	  se += (out(pattern, c) - output(c)) * (out(pattern, c) - output(c));
	}
      }
      se /= in.dim(0);
    }
    return 100 * se / (outputDimension * (outmax - outmin) * (outmax - outmin)) ;
  }

//===========================================================================
/*!
 *  \brief Calculates the mean squared error between the output and the target
 *         vector.
 *
 *  Measures the euklidian distance between the model output \em model(in),
 *  calculated from the input vector \em in, and the target vector \em out. 
 *  The result is then normalized to the number output neurons.
 *  Consider the case of a N-dimensional output vector, i.e. a neural network 
 *  with \em N output neurons, and a set of \em P patterns. In this case the 
 *  function calculates 
 *  \f[
 *      E = \frac{1}{N} \sum_{p=1}^P\sum_{i=1}^N(model(in)_{ip} - 
 *      out_{ip})^{2}
 *  \f]
 *
 *      \param  input Input vector for the model.
 *      \param  target Target vector.
 *      \return The mean squared error \em E.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      C. Igel, M. Toussaint, 2001-09-13:<BR>
 *      Normalized to the number of output neurons instead of
 *      the number of patterns.
 *
 *  \par Status
 *      stable
 *
 */  
  double meanSquaredError(const Array<double> &input, const Array<double> &target) 
  {
    return squaredError(input, target) / target.nelem(); 
  }

//===========================================================================
/*!
 *  \brief Calculates the squared error between the output and the target
 *         vector.
 *
 *  Measures the euklidian distance between the model output \em model(in),
 *  calculated from the input vector \em in, and the target vector \em out. 
 *  Consider the case of a N-dimensional output vector, i.e. a neural network 
 *  with \em N output neurons, and a set of \em P patterns. In this case the 
 *  function calculates 
 *  \f[E = \sum_{p=1}^P\sum_{i=1}^N(model(in)_{ip} - out_{ip})^{2}\f]
 *
 *      \param  input Input vector for the model.
 *      \param  target Target vector.
 *      \return The squared error \em E.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
  double squaredError(const Array<double> &input, const Array<double> &target) {
    double se = 0;
    if(input.ndim() == 1) {
      Array<double> output(target.dim(0));
      model(input, output);
      for(unsigned c = 0; c < target.dim(0); c++) {
	se += (target(c) - output(c)) * (target(c) - output(c));
      }
    } else {
      Array<double> output(target.dim(1));
      for(unsigned pattern = 0; pattern < input.dim(0); ++pattern) {
	model(input[pattern], output);
	for(unsigned c = 0; c < target.dim(1); c++) {
	  se += (target(pattern, c) - output(c)) * (target(pattern, c) - output(c));
	}
      }
    }
    return se;
  }

//===========================================================================
/*!
 *  \brief Caculates the variance of data in each column of a data set.
 *
 *  Assume a data set with \f$P\f$ rows (patterns) and \f$N\f$ 
 *  columns (no. of input neurons), then the variance vector \em v
 *  for \em data is defined by:
 *
 *  \f$
 *  v = ( v_i ),\ \mbox{for\ } i = 1, \dots, N
 *  \f$
 *
 *  \f$
 *  v_i = \frac{1}{P}\sum_{p=1}^P(data_{ip})^{2} -
 *        \left(\frac{1}{P}\sum_{p=1}^Pdata_{ip} \right) ^2
 *  \f$
 *
 *      \param  data Set of data (matrix of input patterns).
 *      \param  v The vector with the calculated variance values.  
 *      \return None.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
  void variance(const Array<double> &data, Array<double> &v) {
    Array< double > m;
    m.resize(data.dim(1));
    v.resize(data.dim(1));
    m = 0.;
    v = 0.;
    for(unsigned i = 0; i < data.dim(0); ++i ) {
      m += data[i];
      v += data[i] * data[i];
    }

    m /= double(data.dim(0));
    v /= double(data.dim(0));
    v -= m * m;
  }    

//===========================================================================
/*!
 *  \brief Caculates the error based on "the winner takes it all" strategy.
 *
 *  This error measure is used for classification problems.
 *  For each pattern only the output neuron with the best (maximum) value
 *  is considered and set to 1, all other output neurons are set
 *  to zero. This result is compared with the target vectors
 *  and then the number of different vectors (= errors) is counted.
 *  Finally, the number of errors is normalized with the number of patterns.
 *
 *      \param  in Input vector for the model.
 *      \param  out Target vector, that must be in binary format.
 *      \return The number of errors, normalized with the number of patterns.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
  double wta(const Array<double> &in, const Array<double> &out) {
    double ce = 0;
    if(in.ndim() == 1) {
      Array<double> output(out.dim(0));
      model(in, output);
      takeAll(output);
      if(out != output) ce++;
    } else {
      Array<double> output(out.dim(1));
      for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	model(in[pattern], output);
	takeAll(output);
	for(unsigned i = 0; i < output.nelem(); i++) {
	  if(out[pattern](i) != output(i)) {
	    ce++;
	    break;
	  }
	}
      }
      ce /= in.dim(0);
    }
    return ce;
  }

//===========================================================================
/*!
 *  \brief Uses thresholds for transforming output and target vectors
 *         into binary vectors and then calculating the error. 
 *
 *  Used for qualification strategies. The target vectors are transformed
 *  into binary vectors by using a threshold of "0.5". All target neurons
 *  greater than "0.5" are set to "1", all other target neurons
 *  are set to zero. The output vectors are transformed into binary
 *  vectors in the same way, but by using the parameter \em t, you
 *  can define a "buffer interval" of size \f$2 \ast t\f$. 
 *  All output neurons with a value
 *  of 0.5 - \em t or less are set to zero, all output neurons with a value
 *  greater than 0.5 + \em t are set to 1, output neurons with a
 *  value in the interval \f$2 \ast t\f$ are ignored. 
 *  This make for these values counted as errors, because normally
 *  the unmodified values of the output neurons are not "0" or "1".
 *  After this transformation of the vectors the number of
 *  vectors is counted for which the used model calculates
 *  another output than given by the corresponding target vector.
 *  Then the number of different vector-pairs is normalized
 *  by the number of input patterns. Given a number of \f$P\f$
 *  patterns the error is calculated by:
 *
 *  \f$
 *    E = \frac{1}{P}\ |A|\ \mbox{with\ }
 *    A = \{ in_i\ |\ binary(model(in_i)) \neq binary(out_i),\ i = 1, \dots, 
 *    P\}
 *  \f$
 *
 *  \em binary is used here for the notation of the transformation
 *  of the original vectors to binary vectors.
 *
 *  \par Example
 *
 *  \f$
 *    \mbox{model(in)\ } = 
 *    \left( 
 *    \begin{array}{c}
 *    0.1\\ 0.9\\ 0.5\\ 0.4\\ 0.3\\ 0.7\\ 0.5\\ 0.1\\ 0.8\\ 0.6\\
 *    \end{array} 
 *    \right) \Rightarrow
 *    \left(
 *    \begin{array}{c}
 *    0\\ 1\\ 0\\ 0\\ 0\\ 1\\ 0\\ 0\\ 1\\ 1\\
 *    \end{array}
 *    \right)
 *    \mbox{\ \ \ \ \ out\ } =
 *    \left( 
 *    \begin{array}{c}
 *    0.6 \\ 0.6 \\ 0.1 \\ 0.2 \\ 0.1 \\ 0.9 \\ 0.7 \\ 0.9 \\ 0.1 \\ 0.8\\    
 *    \end{array} 
 *    \right) \Rightarrow
 *    \left(
 *    \begin{array}{c}
 *    1\\ 1\\ 0\\ 0\\ 0\\ 1\\ 1\\ 1\\ 0\\ 1\\
 *    \end{array}
 *    \right)
 *    \f$
 *
 *    \f$
 *    \mbox{model(in)\ } = 
 *    \left( 
 *    \begin{array}{c}
 *    0.1 \\ 0.9 \\ 0.5 \\ 0.4 \\ 0.3 \\ 0.7 \\ 0.5 \\ 0.1 \\ 0.8 \\ 0.6
 *    \end{array} 
 *    \right) \Rightarrow
 *    \left(
 *    \begin{array}{c}
 *    0 \\ 1 \\ 0.5 \\ 0.4 \\ 0 \\ 0.7 \\ 0.5 \\ 0 \\ 1 \\ 0.6    
 *    \end{array}
 *    \right)
 *    \mbox{\ \ \ \ \ out\ } =
 *    \left( 
 *    \begin{array}{c}
 *    0.6 \\ 0.6 \\ 0.1 \\ 0.2 \\ 0.1 \\ 0.9 \\ 0.7 \\ 0.9 \\ 0.1 \\ 0.8
 *    \end{array} 
 *    \right) \Rightarrow
 *    \left(
 *    \begin{array}{c}
 *    1 \\ 1 \\ 0 \\ 0 \\ 0 \\ 1 \\ 1 \\ 1 \\ 0 \\ 1
 *    \end{array}
 *    \right)\\
 *  \f$
 *
 *  The first example shows the binary transformation for a
 *  value for parameter \em t of "0" resulting in
 *  4 differences between the output and the target vector, 
 *  the second shows the same transformation for a value of "0.2",
 *  resulting in twice as many differences.
 *
 *      \param  in Input vector for the model.
 *      \param  out Target vector.
 *      \param  t Defines an interval of length 2\em t, output values
 *                in this interval are not transformed into
 *                binary values. 
 *      \return 
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
  double binaryCriterion(const Array<double> &in, const Array<double> &out, double t = 0.) {
    double ce = 0;
    
    Array<double > o = out;
    binary(o, 0);

    if(in.ndim() == 1) {
      Array<double> output(out.dim(0));
      model(in, output);
      binary(output, t);
      for(unsigned i = 0; i < output.nelem(); i++) {
	if(o(i) != output(i)) {
	  ce++;
	  break;
	}
      }
    } else {
      Array<double> output(out.dim(1));
      for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	model(in[pattern], output);
	binary(output, t);
	for(unsigned i = 0; i < output.nelem(); i++) {
	  if(o[pattern](i) != output(i)) {
	    ce++;
	    break;
	  }
	}
      }
      ce /= (double) in.dim(0);
    }
    return ce;
  }

//===========================================================================
/*!
 *  \brief Calculates the cross entropy error.
 *
 *  The cross entropy function for \em N patterns and 
 *  \em C>1 class-dimensions within the output vector is calculated via
 *  \f[
 *      E = - \sum_{i=1}^N \sum_{k=1}^C \left\{tar^i_k\cdot \ln \frac{\exp{(model_k(in^i))}}
 *                          {\sum_{k^{\prime}=1}^C \exp{(model_{k^{\prime}}(in^i))}}\right\}
 *  \f] 
 *  respectively for only one single output dimension and two classes via
 *  \f[
 *      E = - \sum_{i=1}^N \left\{tar^i\cdot \ln model(in^i) + (1-tar^i) \ln 
 *          (1-model(in^i))\right\}
 *  \f]
 *
 *      \param  in Input-Matrix for the model.
 *      \param  target Target vector.
 *      \return The error \em E.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *    Revision 2003/06/03 (S. Wiegand): 
 *    softmax activation introduced, bug fixed for multiple class problems
 *
 *  \par Status
 *      stable
 *
 */
  double crossEntropy(const Array<double> &in, const Array<double> &target) {
    const double minLog=-730; // variable is only used to circumvent numerical problems with log(.)
     double ce=0;
    double smaxsum;
    if (outputDimension>1)
      // more then one output neuron
      if(in.ndim() == 1) {
	Array<double> output(target.dim(0));
	model(in, output);
	smaxsum=smaxnorm(output);
	output/=smaxsum;
	for(unsigned c = 0; c < target.dim(0); c++)
	  if(target(c) != 0) 
	    if(output(c) != 0.)
	      ce -= target(c) * log( output(c) ) ;
	    else
	      ce -= minLog;
      } else {
	Array<double> output(target.dim(1));
	for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	  model(in[pattern], output); 
// 	  cout << endl;
// 	  writeArray(output,cout);
	  smaxsum=smaxnorm(output);
// 	  cout << "Maxsum " << smaxsum << endl;
	  output/=smaxsum;
// 	  writeArray(output,cout);
	  for(unsigned c = 0; c < target.dim(1); c++)
	    if(target(pattern, c) != 0)
	      if(output(c) != 0.) 
		ce -= target(pattern, c) * log( output(c) );
	      else
		ce -= minLog;	
	}
      }
    else{
      // only one output neuron
      if(in.ndim() == 1) {
	// one pattern
	Array<double> output(1);
	model(in, output);
	if ((output(0)!=0) && (1. - output(0)!=0))
	  ce -= (target(0)*log(output(0)) + (1. - target(0)) * log(1. - output(0)));
	else if (output(0)!=target(0))
	  ce -= minLog;
      } else {
	// more then one pattern
	Array<double> output(1);
	for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	  model(in[pattern], output);
	  if ((output(0) != 0) && (1. - output(0)!=0))
	    ce -= target(pattern, 0) * log(output(0)) + (1. - target(pattern, 0)) * log(1. - output(0));
	  else if (output(0) != target(pattern, 0))
	    ce -= minLog;	  
	}
      }
    }
    return ce;
  }

 private:

  double smaxnorm(Array<double> &output) {  
    const static double maxexp=1e300;

    double sum = 0;
    if(output.ndim() == 1) {
      for(unsigned i = 0; i < output.nelem(); i++) {
	if(output(i)<301)
	  output(i)=exp(output(i));
	else output(i)=maxexp/output.nelem();
	sum += output(i);
      }
    } else {
      exit(1);
    }
    return sum;
  }

 protected:

  // First searches for the best (maximum) output value,
  // then sets this value to 1 and all others to 0. 
  //
  void takeAll(Array<double> &a) {
    unsigned i, bestIndex = 0;
    double best = a(bestIndex);
    for(i = 1; i < a.nelem(); i++) {
      if(a(i) > best) {
	bestIndex = i;
	best = a(bestIndex);
      }
    }
    for(i = 0; i < a.nelem(); i++) {
      if(i == bestIndex) a(i) = 1.;
      else a(i) = 0.;
    }
  }

  // Used to mark output neurons as set/unset.
  // Normally the threshold for a neuron to be
  // set is 0.5 and higher, but with parameter t
  // a buffer interval between set and unset neurons
  // is defined. All neurons lying in this interval
  // are ignored.
  // 
  void binary(Array<double> &a, double t) {
    if (a.ndim()==1)
      for(unsigned i = 0; i < a.nelem(); i++) {
	if(a(i) > (0.5 + t)) a(i) = 1.;
	if(a(i) <= (0.5 - t)) a(i) = 0.;
      }
    else
      for(unsigned i = 0; i < a.dim(0); i++)
	for (unsigned j=0; j<a.dim(1); j++){
	  if(a(i,j) > (0.5 + t)) a(i,j) = 1.;
	  if(a(i,j) <= (0.5 - t)) a(i,j) = 0.;
	}
  }
};

#endif














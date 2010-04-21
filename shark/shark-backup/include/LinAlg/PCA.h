//===========================================================================
/*!
 *  \file PCA.h
 *
 *  \brief Methods for the "Principal Component Analysis" class.
 *
 *  This file offers methods for the \em Principal \em Component
 *  \em Analysis class, that is used to compress data for a
 *  better visualization and analysis.
 *
 *  \author  M. Kreutz
 *  \date    1998-10-15
 *
 *  \par Copyright (c) 1998-2000:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>
 *
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: PCA.h,v $<BR>
 *      $Revision: 2.3 $<BR>
 *      $Date: 2004/11/19 17:09:19 $
 *
 *  \par Changes:
 *      $Log: PCA.h,v $
 *      Revision 2.3  2004/11/19 17:09:19  shark-admin
 *      documentation added
 *
 *      Revision 2.2  2004/11/19 15:41:33  shark-admin
 *      Methods setState(...) and getState() for PCA objects added
 *
 *      Revision 2.1  2004/06/16 15:22:20  saviapbe
 *
 *      Some bugs in the doxygen documentation were removed.
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.5  2002/02/06 14:28:46  rudi
 *      Error in method "rtransform" removed, new method for extraction of eigenvalues added.
 *
 *      Revision 1.4  2001/11/30 14:11:50  rudi
 *      doxygen comments added.
 *
 *
 *  This file is part of LinAlg. This library is free software;
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

#ifndef PCA_H
#define PCA_H

#include "Array/ArrayIo.h"
#include "Array/Array2D.h"


// The command in the following comment is used for the
// doxygen documentation program. doxygen is forced to include the
// source code of the named example programs into the documentation,
// that then can be linked directly.
// This "dirty trick" is necessary, because to put "\example [filename]"
// directly into method documentations has some unaesthetic side effects.
//
/*! \example pca_test.cpp */



//===========================================================================
/*!
 *  \brief "Principal Component Analysis" class for data compression.
 *
 *  The \em Principal \em Component \em Analysis also known as
 *  \em Karhunen-Loeve-Transformation takes a symmetric
 *  \f$ n \times n \f$ matrix \f$ A \f$ and uses its decomposition 
 *
 *  \f$
 *      A = \Gamma \Lambda \Gamma^T
 *  \f$ 
 * 
 *  where \f$ \Lambda \f$ is the diagonal matrix of eigenvalues
 *  of \f$ A \f$ and \f$ \Gamma \f$ is the orthogonal matrix
 *  with the corresponding eigenvectors as columns. 
 *  \f$ \Lambda \f$ then defines a successive orthogonal rotation, 
 *  that maximizes
 *  the variances of the coordinates, i.e. the coordinate system
 *  is rotated in such a way that the correlation between the new
 *  axes becomes zero. If there are \f$ p \f$ axes, the first
 *  axis is rotated in a way that the points on the new axis
 *  have maximum variance. Then the remaining \f$ p - 1 \f$
 *  axes are rotated such that a another axis covers a maximum
 *  part of the rest variance, that is not covered by the first
 *  axis. After the rotation of \f$ p - 1 \f$ axes, 
 *  the rotation destination of axis no. \f$ p \f$ is fixed.
 *  An application for \em PCA is the reduction of dimensions
 *  by skipping the components with the least corresponding 
 *  eigenvalues/variances.
 *
 *  \example pca_test.cpp
 *
 *  This example (you can find the executable in "$SHARKDIR/LinAlg/examples"
 *  after you have called \f$make\ examples\f$) shows the application of the 
 *  PCA method on a 2-dimensional, randomly created data set. The data
 *  points of the original data and the transformed data points
 *  are saved in four files, named \em pre.dat, \em pca.dat, \em white.dat and
 *  \em red.dat, after the execution of the program.
 *  File \em pre.dat contains the original data, \em pca.dat
 *  the resulting data after performing a PCA, file \em white.dat
 *  the data that come to existence after performing a PCA with \em whitening
 *  (s.a. method #setWhitening) and \em red.dat is the 
 *  original data after performing a PCA and a dimension reduction.
 *  The mean value of all transformed data is moved to the zero point. 
 *  Additionally, the two files \em pca_rt.dat and \em white_rt.dat
 *  contain the data of the files \em pca.dat and \em white.dat
 *  respectively on which a retransformation was performed. The data
 *  in these files should be identical to the original data in
 *  \em pre.dat.
 *  Use the program \em gnuplot to visualize the results. Start
 *  \em gnuplot by typing \f$ gnuplot \f$ and then enter 
 *  \f$ plot\ "pre.dat",\ "pca.dat",\ "white.dat",\ "red.dat" \f$.
 *  The opening window will give you a good overview about
 *  the effect of the different PCA-transformations.
 *  Then enter \f$ plot\ "pre.dat",\ "pca.dat",\ "pca\_rt.dat" \f$
 *  and \f$ plot\ "pre.dat",\ "white.dat",\ "white\_rt.dat" \f$.
 *  This shows you, that the retransformation is working. The
 *  data points of \em pre.dat and \em pca_rt.dat/white_rt.dat
 *  should lie one upon the other.
 *
 *  \author  M. Kreutz
 *  \date    1998-10-15
 *
 *  \warning none
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 *
 */
class PCA
{
  public:

    //! Constructs an new "PCA" object with default values for internal
    //! variables.
    PCA ( );

    //! Destructs an "PCA" object.
    ~PCA ( );

    //! Resets internal variables of a "PCA" object to initial values.
    void reset( );

    //! Updates the transformation matrix.
    void update( );

    //! Calculates the estimates of the mean and the covariance
    //! matrix for the sample data "v".
    void train( const Array< double >& );

    //! Transformation of original data "v" to compressed data "w".
    void transform( const Array< double >&, Array< double >&, unsigned = 0 );

    //! Retransformation of compressed data "v" to original data "w".
    void rtransform( const Array< double >&, Array< double >&, unsigned = 0 );

    //! If set to "true" a matrix will be created in method
    //! #update, that transforms a given sequence of data
    //! to white noise. 
    void setWhitening( bool flagA );

    //! If set to "true" a matrix will be created in method
    //! #update, that will move the data in a way, that
    //! the mean value of the data lies in the zero-point.  
    void setRemoveMean( bool flagA );

    //! Returns the transformation matrix.
    Array< double > transMat( );

    //! Returns the eigenvalues evaluated during the calculation
    //! of the transformation matrix.
    Array< double > eigVal( );

    //! Returns the number of training vectors.
    unsigned getNoTrainVecs() { return countM; }

    //! Returns a pointer to data of transformation matrix.
    //! User must ensure, that "update" has been called before.
    const double*  transMatPtr() const { return transMatM.begin(); }

    //! Returns a pointer to data of eigenvalue matrix.
    //! User must ensure, that "update" has been called before.
    const double*  eigValPtr() const   { return dvecL.begin(); }


    //! returns state of PCA object.
    const void getState(Array2D< double > &cov, Array< double > &m, unsigned &n) 
    {
      cov = covarM;
      m = meanM;
      n = countM;
    }

    //! sets state of PCA object.
    void setState(Array2D< double > &cov, Array< double > &m, unsigned n) 
    {  
      covarM = cov;
      meanM = m;
      countM = n;
      modifiedM = true;
      update();
    }


  private:

    // If set to "true", the transformation matrix must be
    // updated. 
    bool modifiedM;

    // If set to "true" a matrix will be created in method
    // "update()", that transforms a given sequence of data
    // to white noise. See also documentation of 
    // method "setWhitening()".
    bool whiteningM;

    // If set to "true" a matrix will be created in method
    // "update()", that will move the data in a way, that
    // the mean value of the data lies in the zero-point. 
    bool removeMeanM;

    // If data used for training consists of several
    // vectors then a training is performed for each
    // single vector. This flag counts the number of
    // "vector" trainings that were performed.
    unsigned countM;

    // Stores information for the calculation of mean values.
    Array< double > meanM;

    // Stores information for the calculation of the covariance matrix.
    Array2D< double > covarM;

    // The transformation matrix.
    Array2D< double > transMatM;

    // Used to store the eigenvalues when calculating the 
    // transformation matrix.
    Array< double > dvecL;

};

#endif // PCA_H



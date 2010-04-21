//===========================================================================
/*!
 *  \file linearClassifer_test.cpp
 *
 *
 *  \par Copyright (c) 1998-2003:
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
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: linearClassifier_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: linearClassifier_test.cpp,v $
 *      Revision 2.1  2004/03/05 17:36:07  shark-admin
 *      cast bug fixed in  error += std::abs(...)
 *
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
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
#include "Array/ArrayIo.h"
#include "LinAlg/LinearClassifier.h"

using namespace std;

// Number of training vectors:
const unsigned no_train_vecs = 10;

// Number of test vectors:
const unsigned no_test_vecs = 2;

 
// Values for the training vectors:
double train_vec_values[ 2 * no_train_vecs ] = {
    1.79, 51.3,
    1.52, 40.9,
    1.64, 73.0,
    1.67, 64.5,
    1.82, 65.9,
    1.81, 125.3,
    1.62, 100.3,
    1.82, 100.1,
    2.01, 90.3,
    1.8, 51.4
};
 
// Numbers of classes, to which the training vectors belong:
unsigned class_no_values[ no_train_vecs ] = {
    1,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    1
};

// Values for the test vectors:
double test_vec_values[ 2 * no_test_vecs ] = {
    1.80, 65.9,
    1.49, 38.5,
}; 
 
int main()
{
    unsigned          curr_indiv( 0 ),                // Index of current 
                                                      // training vector.
                      class_no;                       // Resulting class number
                                                      // for current test 
                                                      // vector.
    double            error(0.);                      // Cross validation 
                                                      // error.
    Array< double >   train_vecs( no_train_vecs, 2 ), // Training vectors.
                      test_vecs( no_test_vecs, 2 );   // Test vectors.
    Array< unsigned > class_nos( no_train_vecs ),     // The corresponding
                                                      // class numbers.
                      results( no_train_vecs );       // Classification results
                                                      // during cross 
                                                      // validation.
 
 
    // Create a Linear Classifier object:
    LinearClassifier lc;
 
    // Assign values for training and test vectors to storage variables:
    for ( curr_indiv = 0; curr_indiv < no_train_vecs; curr_indiv++ ) 
    {
        train_vecs( curr_indiv, 0 ) = train_vec_values[ curr_indiv * 2 ];
        train_vecs( curr_indiv, 1 ) = train_vec_values[ curr_indiv * 2 + 1 ];  
        class_nos( curr_indiv ) = class_no_values[ curr_indiv ]; 
        if ( curr_indiv < no_test_vecs )
	{
            test_vecs( curr_indiv, 0 ) = test_vec_values[curr_indiv * 2 ];
            test_vecs( curr_indiv, 1 ) = test_vec_values[curr_indiv * 2 + 1 ]; 
        }
    }
 
    // Train classificator with the training vectors:
    cout << "Training classificator..."; 
    lc.train( train_vecs, class_nos );  
    cout << "done!" << "\n\n";     

    // Evaluate the cross validation error:
    lc.leaveOneOut( train_vecs, class_nos, results );
    for ( curr_indiv = 0; curr_indiv < no_train_vecs; curr_indiv++ )
    {
#ifdef _WIN32
        error += abs( (int)class_nos( curr_indiv ) - (int)results( curr_indiv ) );
#else   
        error += std::abs( (int)class_nos( curr_indiv ) - (int)results( curr_indiv ) );
#endif
    }
    error /= no_train_vecs;
    cout << "Cross validation error of classificator = " << error << "\n\n";

    // Classification of test vectors and output of result:
    for ( curr_indiv = 0; curr_indiv < no_test_vecs; curr_indiv++ )
    {
        class_no = lc.classify( test_vecs[ curr_indiv ] );
        cout << "Test vector: " << endl;
        writeArray( test_vecs[ curr_indiv ], cout );
        cout << "Resulting class no.: " << class_no << "\n\n";
    }

    return 0;
}






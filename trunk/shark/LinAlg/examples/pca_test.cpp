//===========================================================================
/*!
 *  \file pca_test.cpp
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
 *      $RCSfile: pca_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: pca_test.cpp,v $
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

#include <fstream>
#include "Rng/GlobalRng.h"
#include "Array/ArrayIo.h"
#include "Array/ArrayOp.h"
#include "LinAlg/PCA.h"

using namespace std;

int main()
{
    unsigned i;                     // Counter variable.
    Array< double > data( 100, 2 ); // Original data.
    Array< double > dataTrans;      // Transformed data.
    Array< double > dataReTrans;    // Retransformed data.    
    PCA             pcaV;           // The PCA instance.
    ofstream        pre,            // Storage file for original data.
                    white,          // Storage file for whitened data.
                    white_rt,       // Storage file for data in "white"
                                    // after retransformation.
                    pca,            // Storage file for data on which
                                    // only a PCA is performed.
                    pca_rt,         // Storage file for data in "pca"
                                    // after retransformation.
                    red;            // Storage file for data
                                    // on which a PCA and dimension reduction 
                                    // is performed.


    // Create random data...
    for ( i = 0; i < data.dim( 0 ); i++ ) 
    {
        double x = Rng::gauss( 2, 2 );
        double y = Rng::gauss( -1, 4 );
        data( i, 0 ) = 2 * x + y;
        data( i, 1 ) = x + 2 * y - .3;
    }
    // ...and save it.
    pre.open( "pre.dat" );
    if ( pre ) writeArray( data, pre );
    pre.close( );

    // Perform only PCA and save result:
    pcaV.reset( );
    pcaV.setRemoveMean( true );
    pcaV.setWhitening( false );
    pcaV.train( data );
    pcaV.transform( data, dataTrans );
    pca.open( "pca.dat" );
    if ( pca ) writeArray( dataTrans, pca );
    pca.close( );
    // Output of used transformation matrix:
    cout << "Transformation matrix used for PCA:" << endl;
    writeArray( pcaV.transMat( ), cout );

    // Retransformation of PCA-data, save result:
    pcaV.rtransform( dataTrans, dataReTrans );
    pca_rt.open( "pca_rt.dat" );
    if ( pca_rt ) writeArray( dataReTrans, pca_rt );
    pca_rt.close( );

    // Perform PCA with whitening and save result:
    pcaV.setRemoveMean( true );
    pcaV.setWhitening( true );
    pcaV.train( data );
    pcaV.transform( data, dataTrans );
    white.open( "white.dat" );
    if ( white ) writeArray( dataTrans, white );
    white.close( );
    // Output of used transformation matrix:
    cout << "Transformation matrix used for PCA with whitening:" << endl;
    writeArray( pcaV.transMat( ), cout );

    // Retransformation of whitened PCA-data, save result:
    pcaV.rtransform( dataTrans, dataReTrans );
    white_rt.open( "white_rt.dat" );
    if ( white_rt ) writeArray( dataReTrans, white_rt );
    white_rt.close( );

    // Perform PCA, reduce dimension and save result:
    pcaV.reset( );
    pcaV.setRemoveMean(true);
    pcaV.train( data );
    dataTrans.resize(1, data.dim( 1 ));
    pcaV.transform(data, dataTrans, 1);
    red.open( "red.dat" );
    if ( red ) {
        for ( i = 0; i < dataTrans.dim( 0 ); i++ ) 
        {
	    // For a better visualization the y-value
	    // '0' is added to the reduced data:
            red << dataTrans(i, 0) << " " << 0 << endl;
        }
    }
    red.close( );
    // Output of used transformation matrix:
    cout << "Transformation matrix used for PCA with dimension reduction:" 
         << endl;
    writeArray( pcaV.transMat( ), cout );


    return 0;
}

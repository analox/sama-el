//===========================================================================
/*!
 *  \file rbfn-example.cpp
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
 *      Mixture
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: rbfn-example.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: rbfn-example.cpp,v $
 *      Revision 2.1  2004/03/04 13:45:23  shark-admin
 *      include file iomanip.h exhchanged for iomanip
 *
 *      Revision 2.0  2003/11/28 16:23:11  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of Mixture. This library is free software;
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



#include <Array/ArrayIo.h>
#include <TestData/TimeSeries/DiscreteMackeyGlass.h>
#include <TestData/TimeSeries/IOSamples.h>
#include <ReClaM/BFGS.h>
#include <ReClaM/Rprop.h>
#include <Mixture/RBFN.h>
#include <ReClaM/SquaredError.h>
#include <ReClaM/ErrorMeasures.h>
#include <Rng/RNG.h>
#include <iomanip>

//===========================================================================

class RBFN_Interface : public RBFN, virtual public ModelInterface
{
  public:
    RBFN_Interface( unsigned numInput, unsigned numOutput, unsigned numHidden )
      : RBFN( numInput, numOutput, numHidden )
    {
        ModelInterface::w.resize( b.nelem( ) +
				  A.nelem( ) +
				  m.nelem( ) +
				  v.nelem( ) );
        ModelInterface::dmdw.resize( odim( ),
				   b.nelem( ) +
				   A.nelem( ) +
				   m.nelem( ) +
				   v.nelem( ) );
        ModelInterface::dedw.resize( b.nelem( ) +
				   A.nelem( ) +
				   m.nelem( ) +
				   v.nelem( ) );
        getParams( ModelInterface::w );
    }

    void initialize( const Array< double >& in,
		     const Array< double >& out )
    {
        RBFN::initialize( in, out );
	getParams( ModelInterface::w );
    }

    void model( const Array< double >& in, Array< double >& out )
    {
        setParams( ModelInterface::w );
        recall( in, out );
    }

    void dmodel( const Array< double >& in )
    {
        setParams( ModelInterface::w );
        gradientOut( in, ModelInterface::dedw );
    }
};

class RBFN_MSE : public RBFN_Interface
{
  public:
    RBFN_MSE( unsigned numInput, unsigned numOutput, unsigned numHidden )
      : RBFN_Interface( numInput, numOutput, numHidden ) { }

    double error( const Array< double >& in,
		  const Array< double >& out )
    {
        setParams( ModelInterface::w );
        return mse( in, out );
    }

    double derror( const Array< double >& in,
		 const Array< double >& out,
		 bool returnError = true)
    {
        setParams( ModelInterface::w );
	double e = gradientMSE( in, out, ModelInterface::dedw );
	if (!returnError)
	  e = -1;
	return e;
    }
};

class RBFN_RPROP : public RBFN_MSE,
		   public IRpropPlus,
		   public ErrorMeasures
{
  public:
    RBFN_RPROP( unsigned numInput, unsigned numOutput, unsigned numHidden )
      : RBFN_MSE( numInput, numOutput, numHidden ) { }
};


class RBFN_RPROP2 : public RBFN_MSE,
		    public RpropMinus,
		    public ErrorMeasures
{
  public:
    RBFN_RPROP2( unsigned numInput, unsigned numOutput, unsigned numHidden )
      : RBFN_MSE( numInput, numOutput, numHidden ) { }
};

//===========================================================================

int main( int argc, char **argv )
{
    const unsigned EmbedDim = 5;
    const unsigned TimeLag  = 1;
    const unsigned Horizon  = 1;
    const unsigned NumSkip  = 500;
    const unsigned NumTrain = 1000;
    const unsigned NumTest  = 1000;
    const unsigned NumIter  = 1000;

    unsigned i, t;
    Array< double > inTrain ( NumTrain, EmbedDim );
    Array< double > outTrain( NumTrain, 1        );
    Array< double > inTest  ( NumTest,  EmbedDim );
    Array< double > outTest ( NumTest,  1        );

    DiscreteMackeyGlass mg;
    IOSamples< double > iosamples( mg, EmbedDim, TimeLag, 1, Horizon );

    //
    // create samples for training and test
    //
    for( i = 0; i < NumSkip; ++i ) {
        ArrayReference< double > in;
	ArrayReference< double > out;
	in .copyReference( inTrain [ 0 ] );
	out.copyReference( outTrain[ 0 ] );
	iosamples( static_cast< Array< double >& >( in  ),
		   static_cast< Array< double >& >( out ) );
    }
    for( i = 0; i < NumTrain; ++i ) {
        ArrayReference< double > in;
	ArrayReference< double > out;
	in .copyReference( inTrain [ i ] );
	out.copyReference( outTrain[ i ] );
	iosamples( static_cast< Array< double >& >( in  ),
		   static_cast< Array< double >& >( out ) );
    }
    for( i = 0; i < NumTest; ++i ) {
        ArrayReference< double > in;
	ArrayReference< double > out;
	in .copyReference( inTest [ i ] );
	out.copyReference( outTest[ i ] );
	iosamples( static_cast< Array< double >& >( in  ),
		   static_cast< Array< double >& >( out ) );
    }

    RNG::globalRng.seed( 1 );
    RBFN_RPROP rbfn1( EmbedDim, 1, 10 );
    rbfn1.initialize( inTrain, outTrain );
    rbfn1.initRprop( argc > 1 ? atof( argv[ 1 ] ) : 0.0125 );
    RNG::globalRng.seed( 1 );
    RBFN_RPROP2 rbfn2( EmbedDim, 1, 10 );
    rbfn2.initialize( inTrain, outTrain );
    rbfn2.initRprop( argc > 1 ? atof( argv[ 1 ] ) : 0.0125 );

    for( t = 0; t < NumIter; ++t ) {
      //rbfn1.train( inTrain, outTrain );
        rbfn1.rprop( inTrain, outTrain );
	std::cout << t << '\t'
	     << rbfn1.mse( inTrain, outTrain ) << '\t'
	     << rbfn1.mse( inTest,  outTest  ) << '\t';
        rbfn2.rprop( inTrain, outTrain );
	std::cout << rbfn2.mse( inTrain, outTrain ) << '\t'
	     << rbfn2.mse( inTest,  outTest  ) << std::endl;
    }

    return 0;
}

//===========================================================================

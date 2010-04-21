#include "fifo.h"

FIFO::FIFO() {

}

FIFO::~FIFO() {

}

FIFO::FIFO(unsigned maxSize, unsigned inputDim, unsigned targetDim) {
	maximumSize = maxSize;
	inputDimension = inputDim;
	targetDimension = targetDim;
	currentSize = 0;

	input.resize( maxSize, inputDimension, false );
	target.resize( maxSize, targetDimension, false );

	initArrayToZero( input );
	initArrayToZero( target );
}

void FIFO::append(vector<double>& inputData, vector<double>& targetData) {
	//assert( inputData.size() == inputDimension 
	//	&& targetData.size() == targetDimension );

	unsigned i, j;
	
	if( !isFull() ) {
		for( i=0; i<inputDimension; i++ ) {
			input( currentSize, i ) = inputData[ i ];
		}	
		for( i=0; i<targetDimension; i++ ) {
			target( currentSize, i ) = targetData[ i ]; 
		}
		
		currentSize++;
	}
	else {
		for( i=0; i<maximumSize-1; i++ ) {
			for( j=0; j<inputDimension; j++ ) {
				input( i, j ) = input( i+1, j );
			}
			for( j=0; j<targetDimension; j++ ) {
				target( i, j ) = target( i+1, j );
			}
		}

		for( i=0; i<inputDimension; i++ ) {
			input( maximumSize-1, i ) = inputData[ i ]; 
		}
		for( i=0; i<targetDimension; i++ ) {
			target( maximumSize-1, i ) = targetData[ i ]; 
		}
	
	}
}

void FIFO::append( Array<double>& inputData, Array<double>& targetData) {
	//assert( inputData.dim(0)==targetData.dim(0)
	//	&& inputData.dim(1)==inputDimension
	//	&& targetData.dim(1)==targetDimension );

	unsigned i, j;
	unsigned nInput = inputData.dim(0);

	for( i=0; i<nInput; i++ ) {
		vector< double > tempInput(inputDimension);
		vector< double > tempTarget(targetDimension);

		for( j=0; j<inputDimension; j++ ) {
			tempInput[ j ] = inputData( i, j );
		} 
		for( j=0; j<targetDimension; j++ ) {
			tempTarget[ j ] = targetData( i, j );
		}

		append( tempInput, tempTarget );
	}
}

void FIFO::replace(vector<double>& inputData, vector<double>& targetData, unsigned index) {
	//assert( inputData.size() == inputDimension 
	//	&& targetData.size() == targetDimension 
	//	&& index < maximumSize );

	unsigned i;
	for( i=0; i<inputDimension; i++ ) {
		input( index, i ) = inputData[ i ];
	}
	for( i=0; i<targetDimension; i++ ) {
		target( index, i ) = targetData[ i ];
	}
}

unsigned FIFO::getCurrentSize() {
	return currentSize;
}

unsigned FIFO::getMaximumSize() {
	return maximumSize;
}

void FIFO::setCurrentSize( unsigned currSize ) {
	currentSize = currSize;
}

void FIFO::setMaximumSize( unsigned maxSize ) {
	maximumSize = maxSize;
}

bool FIFO::isFull() {
	return ( currentSize == maximumSize? true: false );
}

Array<double> FIFO::getInput() {
	return input;
}

Array<double> FIFO::getInputAt(int index) {
	Array<double> arr(1, inputDimension);
	for(int i=0; i<inputDimension; i++) {
		arr(0,i) = input(index, i);
	}
	return arr;
}

Array<double> FIFO::getTargetAt(int index) {
        Array<double> arr(1, targetDimension);
        for(int i=0; i<targetDimension; i++) {
                arr(0,i) = target(index, i);
        }
        return arr;
}


Array<double> FIFO::getTarget() {
	return target;
}

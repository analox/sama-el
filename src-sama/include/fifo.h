#ifndef _FIFO_
#define _FIFO_
#include <vector.h>
#include "Array/Array.h"

class FIFO {
	private:
		Array< double > input;
		Array< double > target;
		unsigned inputDimension;
		unsigned targetDimension;
		unsigned currentSize;
		unsigned maximumSize;
		
	public:
		FIFO();
		~FIFO();
		FIFO(unsigned maxSize, unsigned inputDim, unsigned targetDim);
		void append(vector<double>& inputData, vector<double>& targetData);
		void append(Array<double>& inputData, Array<double>& targetData);
		void replace(vector<double>& inputData, vector<double>& targetData, unsigned index);
		Array<double> getInput();
		Array<double> getTarget();
		unsigned getCurrentSize();
		unsigned getMaximumSize();
		void setCurrentSize( unsigned currSize );
		void setMaximumSize( unsigned maxSize );
		bool isFull();		

		Array<double> getInputAt(int i);
		Array<double> getTargetAt(int i);

                template <class Type>
                static void initArrayToZero( Array< Type >& arr ) {
                        unsigned i, j;
                        for ( i = 0; i < arr.dim( 0 ); i++ ) {
                                for ( j = 0; j < arr.dim( 1 ); j++ ) {
                                        arr( i, j ) = 0;
                                }
                        }
                }

};


#endif

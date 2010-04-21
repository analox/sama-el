#include <iostream.h>
#include "fifo.h"

int main(int argc, char** argv) {
	unsigned i, j;
	unsigned inDim = 2, outDim = 2, nInput = 20;
	FIFO myList( nInput, inDim, outDim );
	
	Array< double > newInput(3, inDim);
	Array< double > newTarget(3, outDim);

	for(i=0; i<3; i++) {
		for(j=0; j<inDim; j++) {
			newInput(i,j) = 2;	
		}
		for(j=0; j<outDim; j++) {
			newTarget(i,j) = i;	
		}
	}

	vector<double> newRepl(inDim);
	vector<double> newReplFit(outDim);

	for(j=0; j<inDim; j++) {
		newRepl[j] = 2;	
	}
	for(j=0; j<outDim; j++) {
		newReplFit[j] = 5;	
	}

	myList.replace( newRepl, newReplFit, 5 );

	myList.append(newInput, newTarget);

	for(i=0; i<nInput; i++) {
		for(j=0; j<inDim; j++) cout << (myList.getInput())( i, j) << "\t";
		for(j=0; j<outDim; j++) cout << (myList.getTarget())( i, j) << "\t";

		cout<<endl;
	}

	for(j=0; j<inDim; j++) cout << (myList.getTargetAt(5))( 0, j) << "\t";
	cout << endl;

	return 0;
}

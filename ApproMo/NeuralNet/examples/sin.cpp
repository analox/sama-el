#include "Population.h"
#include "Array.h"
#include "EvoNet.h"
int main( int argc, char **argv )
{
  unsigned inDimension = 2;
  unsigned outDimension = 2;
  unsigned numberOfHiddenNeurons = 15;
  unsigned numberOfNeurons = numberOfHiddenNeurons+inDimension+outDimension;
  
  Array< int > cmat( numberOfNeurons, numberOfNeurons+1 );
  Array<double> wmat( numberOfNeurons, numberOfNeurons+1 );
  //Initialize a feed-forward network
  
  for ( unsigned ii = inDimension; ii < inDimension + numberOfHiddenNeurons; ii++ ){
    for ( unsigned jj = 0; jj < inDimension ; jj++ ){
      cmat( ii, jj ) = 1;
      wmat( ii, jj) = Rng::uni(-0.01,0.01);
    }
  }
  // Set bias values:
  for ( unsigned ii = inDimension; ii < inDimension + numberOfHiddenNeurons+outDimension; ii++ ){
    cmat( ii, cmat.dim( 1 )-1 ) = 1;
    wmat( ii, cmat.dim( 1 )-1) = Rng::uni(-0.01,0.01);;
  }
  
  //Set connections from hidden to output
  
  for ( unsigned ii = inDimension+numberOfHiddenNeurons; ii < inDimension + numberOfHiddenNeurons+outDimension; ii++ ){
    for (unsigned jj =inDimension; jj < inDimension+numberOfHiddenNeurons; jj++ ){
      cmat( ii, jj ) = 1;
      wmat( ii, jj) = Rng::uni(-0.01,0.01);
    }
  }
  
  Database myData;
  EvoNet Network(inDimension, outDimension, cmat, wmat);
  myData.SetMaxArchiveLength(100);
  myData.SetMaxTrainingData(100);
  Network.SetLearnSteps(500);
  MyNet network(inDimension, outDimension, cmat, wmat);
 
  double x=0;
  double y = 0.;
  Array<double> SampleIn(1,inDimension), SampleOut(1,outDimension);
  Array<double> SamIn(100,inDimension), SamOut(100,outDimension);
  
  for(unsigned i=0; i<100; i++){
    SampleIn(0,0) = x;
    SampleIn(0,1) = y;
    SamIn(i,0) = x;
    SamIn(i,1) = y;
    SampleOut(0,0) = sin(2*3.1415926*(x+y));
    SampleOut(0,1) = cos(2*3.1415926*(x+y));
    SamOut(i,0) = sin(2*3.1415926*(x+y));
    SamOut(i,1) = cos(2*3.1415926*(x+y));
    x += 0.01;
    y += 0.02;
    myData.AddTrainingData(SampleIn,SampleOut,i);
  }
  
  printf("N1: Before training: MSE=%f \n", Network.MSE(myData));
  Network.Train_Rprop(myData,1.2,0.5,50,0,0.001);
  printf("N1: After training: MSE=%f \n", Network.MSE(myData));
  
  printf("N2: Before training: MSE=%f \n", network.meanSquaredError(SamIn,SamOut));
  network.initRprop(0.001);
  for(int i=0; i<500; i++)
    network.rprop(SamIn,SamOut, 1.2,0.5,50,0);
  printf("N2: After training: MSE=%f \n", network.meanSquaredError(SamIn,SamOut));
  return 0;
  
}

  
    
  

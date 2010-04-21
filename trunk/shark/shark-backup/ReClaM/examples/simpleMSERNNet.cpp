#include<ReClaM/MSERNNet.h>
#include<ReClaM/Rprop.h>

using namespace std;
#pragma warning(disable: 4250)
//===========================================================================
// some parameters

unsigned iterations=1000;
char* datafile="timeseries";
unsigned episode=1000;
int forecast=5;
bool checkGradient=false;       //if this is true, the gradient calculation is checked
double eps=1e-7;
int connectionCMatrix [3*10*10]={  //this C-Array is later casted to the proper 3-dim array
  1,0,0,0,0,0,0,0,0,0,
  1,1,0,0,0,0,0,0,0,0,
  1,1,1,0,0,0,0,0,0,0,
  1,1,1,1,0,0,0,0,0,0,
  1,1,1,1,1,0,0,0,0,0,
  1,1,1,1,1,1,0,0,0,0,
  1,1,1,1,1,1,1,0,0,0,
  1,1,1,1,1,1,1,1,0,0,
  1,1,1,1,1,1,1,1,1,0,
  1,1,1,1,1,1,1,1,1,1,

  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,

  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,1,1
};


//===========================================================================
// define the regression model

class myNet:public MSERNNet,public RpropPlus{
public:
  myNet():MSERNNet(){
    unsigned dims[3]={3,10,10};
    ArrayReference<int> connectionMatrix(dims,connectionCMatrix,3,3*10*10);
    setStructure(connectionMatrix);
    initWeights(-0.5,0.5);
    initRprop();
  };

  Array<double> &Dedw(){ return dedw; } 
};


//===========================================================================
// learn the Lorentz time series

int main(int argc, char *argv[]){

  unsigned long t=0,i;
  ifstream infile;
  ofstream outfile;

  //load the data
  Array<double> data;
  infile.open(datafile);
  readArray(data,infile);
  infile.close();
  data/=100.; data+=.5; //scale between [0,1] because all neurons of a RNN are sigmoid!

  //define the training and evaluation arrays
  unsigned dims[2]={episode,1};
  ArrayReference<double> trainIn    (dims,&data.elem(0)            ,2,episode);
  ArrayReference<double> trainTarget(dims,&data.elem(0)   +forecast,2,episode);
  ArrayReference<double> evalIn     (dims,&data.elem(1000)         ,2,episode);
  ArrayReference<double> evalTarget (dims,&data.elem(1000)+forecast,2,episode);
  Array<double> trainOut(episode,1),evalOut(episode,1);

  //output the evaluation input/target arrays -- to compare them to evalOut
  outfile.open("input");
  writeArray(evalIn,outfile);
  outfile.close();
  outfile.open("target");
  writeArray(evalTarget,outfile);
  outfile.close();

  //create the network
  myNet net;
  Array<double> Dedw;
  double z;

  while(t++<iterations){
    cout <<t <<"\t";

    //WARNING!!: for a unique error, warmupLength must be >0
    // -- only then, the neuron states are reset to zero
    // in the model routine
    net.includeWarmUp(100);

    //check gradient on the training data
    if(checkGradient){
      net.derror(trainIn,trainTarget);
      Dedw=net.Dedw();
      cout <<"derror - gradient:" <<Dedw;
      net.ModelInterface::derror(trainIn,trainTarget);
      cout <<"deltagrad - gradient:" <<net.Dedw();
      for(z=0,i=0;i<Dedw.nelem();i++){
	Dedw(i)=net.Dedw()(i)-Dedw(i); z+=Dedw(i)*Dedw(i); }
      cout <<"difference:" <<Dedw;
      cout <<"square-sum of difference: " <<z <<"\n";
    }

    //learn on the training data
    net.rprop(trainIn,trainTarget);

    //evaluate on the evaluation data (compare files `output' `input' `target')
    cout <<net.error(evalIn,evalTarget) << endl;
    if(!(t%100)){
      net.model(evalIn,evalOut);
      outfile.open("output");
      writeArray(evalOut,outfile);
      outfile.close();
    }
  }
  return 0;
}

//#pragma implementation

#include "FFOps_OL_neu.h"

void FFOps_OL_neu::initConnections(Individual& Ind, double prob1,
			       double prob2, double prob3){
  unsigned nIn = getNIn(Ind);
  unsigned nHid = getNHid(Ind);
  unsigned nOut = getNOut(Ind);
  Array<int> c = FFManip::initFFOL(nIn, nHid, nOut, prob1, prob2, prob3);
  setStructure(Ind, c, 0);
}

void FFOps_OL_neu::addConnection(Individual& Ind, double low, double high){
  Array<int> c = getConnections(Ind, 0);
  unsigned nIn  = getNIn(Ind);
  unsigned nHid = getNHid(Ind);  

  // bestimme noch freie Verbindungen zwischen Eingangs- und 
  // versteckter Schicht
  Array<unsigned> help(nIn*nHid,2u);
  unsigned max = 0;
  for (unsigned z=nIn; z<nIn+nHid; z++){
    for (unsigned s=0; s<nIn; s++)
      if (c(z,s)==0){
        help(max, 0) = z;
        help(max, 1) = s;
        max++;
      }
  }

  // eine Verbindung hinzufuegen
  if (max > 0){
    Array<double> w = getWeights(Ind, 0);
    unsigned rand = Rng::discrete(0,max-1);
    c(help(rand,0),help(rand,1)) = 1;
    w(help(rand,0),help(rand,1)) = Rng::uni(low, high);
    
    setStructure(Ind,c,0);
    setStructure(Ind,w,0);
  }
}

void FFOps_OL_neu::deleteConnection(Individual& Ind, const double sigma){
  Array<int> c = getConnections(Ind, 0);
  unsigned nIn  = getNIn(Ind);
  unsigned nHid = getNHid(Ind);  

  // bestimme alle besetzen Verbindungen
  // zwischen Eingangs- und versteckter Schicht
  Array<unsigned> help(nIn*nHid,2u);
  unsigned max = 0;
  for (unsigned z=nIn; z<nIn+nHid; z++){
    for (unsigned s=0; s<nIn; s++)
      if (c(z,s)==1){
        help(max, 0) = z;
        help(max, 1) = s;
        max++;
      }
  }

  // eine Verbindung entfernen
  if (max > 0){
    Array<double> w = getWeights(Ind, 0);
    unsigned pos;
    if (sigma<=0){
      pos = 0;
      double min=fabs(w(help(0,0),help(0,1)));
      for (unsigned i=1; i<max; i++)
	if (fabs(w(help(i,0),help(i,1)))<min){
	  min = fabs(w(help(i,0),help(i,1)));
	  pos = i;
	}
    }
    else{
      Array<double> prob(max);  
      for (unsigned i=0; i<max; i++)
      prob(i) = 1.0/(1.0+(w(help(i,0),help(i,1))*w(help(i,0),help(i,1))/(sigma*sigma)));
      //      prob(i) = exp(-w(help(i,0),help(i,1))*w(help(i,0),help(i,1))/sigma/sigma);
      pos = wheelOfFortune(prob);
    }
    c(help(pos,0),help(pos,1)) = 0;
    w(help(pos,0),help(pos,1)) = 0;
    
    setStructure(Ind,c,0);
    setStructure(Ind,w,0);
  }
}

void FFOps_OL_neu::deleteNeuron(Individual& Ind){
  unsigned nHid = getNHid(Ind);
  if (nHid>1){
    unsigned pos = Rng::discrete(0,nHid-1);
    delNeuron(Ind, pos);
  }
}

void FFOps_OL_neu::deleteNeuronGain(Individual& Ind, double sigma){
  unsigned nHid = getNHid(Ind);
  if (nHid>1){
    unsigned nIn = getNIn(Ind);
    Array<double> gain = getGain(Ind);
    unsigned pos;
    if(sigma<=0){
      pos = 0;
      double min=gain(nIn);
      for (unsigned i=1; i<nHid; i++) 
	if (gain(nIn+i) < min){
	  min = gain(nIn+i);
	  pos = i;
	}
    }
    else{
      Array<double> prob(nHid);
      for (unsigned i=0; i<nHid; i++)
	prob(i) = 1.0/(1.0+gain(i+nIn)*gain(i+nIn)/(sigma*sigma));
      
      pos = wheelOfFortune(prob);
    }
    delNeuron(Ind, pos);
  }
}

void FFOps_OL_neu::deleteNeuronActivation(Individual& Ind, double sigma){
  unsigned nHid = getNHid(Ind);
  if (nHid>1){
    unsigned nIn = getNIn(Ind);
    Array<double> activation = getActivation(Ind);
    unsigned pos;
    if(sigma<=0){
      pos = 0;
      double min=activation(nIn);
      for (unsigned i=1; i<nHid; i++) 
	if (activation(nIn+i) < min){
	  min = activation(nIn+i);
	  pos = i;
	}
    }
    else{
      Array<double> prob(nHid);
      for (unsigned i=0; i<nHid; i++)
	prob(i) = 1.0/(1.0+activation(i+nIn)*activation(i+nIn)/(sigma*sigma));
      
      pos = wheelOfFortune(prob);
    }
    delNeuron(Ind, pos);
  }
}


void FFOps_OL_neu::addNeuron(Individual& Ind, double low, double high){
  unsigned nIn = getNIn(Ind);
  unsigned nHid = getNHid(Ind);
  unsigned pos = Rng::discrete(0,nHid);
  insNeuron(Ind, pos);

  // connect new neuron
  Array<int>    c = getConnections(Ind,0);
  Array<double> w = getWeights(Ind, 0);

  // eine Verbindung zwischen Eingabeschicht und neuem verstecktem Neuron
  unsigned vonPos = Rng::discrete(0,nIn-1);
  c(nIn+pos, vonPos) = 1;
  w(nIn+pos, vonPos) = Rng::uni(low, high);
  // Verbindung des Bias
  c(nIn+pos, nIn+pos) = 1;
  w(nIn+pos, nIn+pos) = Rng::uni(low, high);  
  // alle Verbindungen zur Ausgabeschicht
  for (unsigned i=nIn+nHid+1; i<c.dim(0); i++){
    c(i, nIn+pos) = 1;
    w(i, nIn+pos) = Rng::uni(low, high);
  }
  setStructure(Ind, c, 0);
  setStructure(Ind, w, 0);
}

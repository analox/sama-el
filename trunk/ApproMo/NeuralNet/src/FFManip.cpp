//#pragma implementation

#include "FFManip.h"

void FFManip::connectNew(Array<int>& c, Array<double>&w, unsigned pos, double low, double high){
  unsigned vonPos = Rng::discrete(0,pos-1);
  c(pos, vonPos) = 1;
  w(pos, vonPos) = Rng::uni(low, high);
  unsigned nachPos = Rng::discrete(pos+1,c.dim(0)-1);
  c(nachPos, pos) = 1;
  w(nachPos, pos) = Rng::uni(low, high);
}

void FFManip::pullThrough(Array<int> &cOld, Array<int>& c, Array<double>&w, unsigned pos, double low, double high){
  // normale Verbindungen durchschleiden.
  for (unsigned in = 0; in < pos; in++)
    if (cOld(pos,in) == 1)
      for (unsigned out = pos+1; out < cOld.dim(0); out++)
	if (cOld(out, pos) == 1)
	  if (c(out-1,in)==0){
	    c(out-1,in) = 1;
	    w(out-1,in) = Rng::uni(low, high);
	  }

  // bias durchschleifen, falls vorhanden
  if (cOld(pos, pos)==1)
    for (unsigned out = pos+1; out < cOld.dim(0); out++) 
      if ((cOld(out, pos) == 1) & (c(out-1,out-1)==0)){
	c(out-1,out-1)=1;
	w(out-1,out-1)=Rng::uni(low,high);
      }
}

Array<int> FFManip::initFFBulk(unsigned nIn, unsigned nHid, unsigned nOut, double prob){
  Array<int> c(nIn+nHid+nOut,nIn+nHid+nOut);
  c = 0;
  for (unsigned z=nIn; z<c.dim(0); z++)
    for (unsigned s=0; s<=z ; s++)
      if ((s<nIn+nHid) | (s==z))
	if (Rng::uni(0,1)<prob)
	  c(z,s)=1;
  return c;
}

Array<int> FFManip::initFFOL(unsigned nIn, unsigned nHid, unsigned nOut, double prob1, double prob2, double prob3){
  Array<int> c(nIn+nHid+nOut,nIn+nHid+nOut);
  c = 0;
  // block 1
  for (unsigned z=nIn; z<nIn+nHid; z++)
    for (unsigned s=0; s<nIn ; s++)
      if (Rng::uni(0,1)<prob1)
	c(z,s)=1;
  // block 2
  for (unsigned z=nIn+nHid; z<c.dim(0); z++)
    for (unsigned s=nIn; s<nIn+nHid ; s++)
      if (Rng::uni(0,1)<prob2)
	c(z,s)=1;
  // bias
  for (unsigned z=nIn; z<c.dim(0); z++)
    if (Rng::uni(0,1)<prob3)
      c(z,z)=1;

  return c;
}

void FFManip::repairBulk(Array<int>& c, const unsigned nIn, const unsigned nHid, const unsigned nOut){
   // vorwaerts reparieren
  for (unsigned neuron=0; neuron<nIn+nHid; neuron++)
    if (nOutput(c,neuron)==0){
      unsigned inMin = c.dim(1);
      unsigned posMin = 0;
      for (unsigned i=max(nIn, neuron+1); i<nIn+nHid+nOut; i++)
        if (nInput(c,i)<inMin){
          inMin  = nInput(c,i);
          posMin = i;
        }
      c(posMin, neuron) = 1;
    }

  // rueckwaerts reparieren
  for (unsigned neuron=nIn; neuron < nIn+nHid+nOut; neuron++)
    if (nInput(c,neuron)==0){
      unsigned pos = Rng::discrete(0,neuron-1);
      c(neuron, pos) = 1;
    }
}

void FFManip::repairOL(Array<int>& c, const unsigned nIn, 
		       const unsigned nHid, const unsigned nOut){
  // links form input-layer to hidden layer
  for (unsigned neuron = 0; neuron < nIn; neuron++){
    unsigned outlink=0;
    for (unsigned out=nIn; out<nIn+nHid; out++)
      if (c(out, neuron)!=0)
	outlink++;
    if (outlink==0)
      c(Rng::discrete(nIn,nIn+nHid-1), neuron)=1;
  }
  
  for (unsigned neuron = nIn; neuron < nIn+nHid; neuron++){
    unsigned inlink=0;
    for (unsigned in=0; in<nIn; in++)
      if (c(neuron, in)!=0)
	inlink++;
    if (inlink==0)
      c(neuron,Rng::discrete(0,nIn-1))=1;
  }
  
  // links form hidden layer to output-layer
  for (unsigned neuron = nIn; neuron < nIn+nHid; neuron++){
    unsigned outlink=0;
    for (unsigned out=nIn+nHid; out<nIn+nHid+nOut; out++)
      if (c(out, neuron)!=0)
	outlink++;
    if (outlink==0)
      c(Rng::discrete(nIn+nHid,nIn+nHid+nOut-1), neuron)=1;
  }
  
  for (unsigned neuron = nIn+nHid; neuron < nIn+nHid+nOut; neuron++){
    unsigned inlink=0;
    for (unsigned in=nIn; in<nIn+nHid; in++)
      if (c(neuron, in)!=0)
	inlink++;
    if (inlink==0)
      c(neuron,Rng::discrete(nIn,nIn+nHid-1))=1;
  }
}

unsigned FFManip::max(unsigned a, unsigned b){
  if (a>b) return a;
  else return b;
}

unsigned FFManip::min(unsigned a, unsigned b){
  if (a<b) return a;
  else return b;
}

unsigned FFManip::nInput(Array<int>& c, const unsigned N){
  unsigned count = 0;
  for (unsigned s=0; s<N; s++)
    if (c(N,s) == 1) 
      count++;
  return count;
}

unsigned FFManip::nOutput(Array<int>& c, const unsigned N){
  unsigned count = 0;
  for (unsigned z=N+1; z<c.dim(0); z++)
    if (c(z,N) == 1)
      count++;
  return count;
}

Array<unsigned> FFManip::dispensableFFConnection(Array<int>& c, unsigned nIn, unsigned nHid){
  unsigned count =0;
  Array<unsigned> help(c.dim(0)*c.dim(0),2u);
  for (unsigned z=nIn; z<c.dim(0); z++){
    for (unsigned s=0; s<min(z, nIn+nHid); s++)
      if (c(z,s)!=0)
	if ( ( nInput(c,z) > 1 ) && ( ( nOutput(c,s) > 1) || ( s < nIn ) ) ){
	  help(count, 0) = z;
	  help(count, 1) = s;
	  count++;
	}
    if (c(z,z)!=0){
      help(count, 0) = z;
      help(count, 1) = z;
      count++;
    }
  }

  Array<unsigned> help2(count, 2u);
  for (unsigned i=0; i<help2.dim(0); i++){
    help2(i,0) = help(i,0);
    help2(i,1) = help(i,1);
  }

  return help2;
}

Array<unsigned> FFManip::feasibleFFConnection(Array<int>& c, unsigned nIn, unsigned nHid){
  unsigned count =0;
  Array<unsigned> help(c.dim(0)*c.dim(0),2u);
  for (unsigned z=nIn; z<c.dim(0); z++){
    for (unsigned s=0; s<min(z,nIn+nHid); s++)
      if (c(z,s)==0){
	help(count, 0) = z;
	help(count, 1) = s;
	count++;
      }
    if (c(z,z)==0){
      help(count, 0) = z;
      help(count, 1) = z;
      count++;
    }
  }

  Array<unsigned> help2(count, 2u);
  for (unsigned i=0; i<help2.dim(0); i++){
    help2(i,0) = help(i,0);
    help2(i,1) = help(i,1);
  }

  return help2;
}

//#pragma implementation

#include "FFOps_Basic_neu.h"

Individual FFOps_Basic_neu::createEmpty(const unsigned in, const unsigned hid, 
			 const unsigned out){
  Individual Ind=Ops_Basic::createEmpty(in, hid, out, 1);
  ChromosomeT<double> w;
  Ind.append(w);
  //  Ind.append(w);
  return Ind;
  }

// void FFOps_Basic_neu::setGain(Individual& Ind, Array<double> gain){
//   if(gain.dim(0) != getNNeurons(Ind)){
//     cerr << "error setGain: dimension of gain vector (" << gain.dim(0) 
// 	 << ") and number of neurons (" << getNNeurons(Ind) << ") do not coincide." << endl;
//     exit(EXIT_FAILURE);
//   }

//   ((ChromosomeT<double>&) (Ind[3])).resize(getNNeurons(Ind));
//   for (unsigned i=0; i<gain.dim(0); i++)
//     ((ChromosomeT<double>&) (Ind[3]))[i] = gain(i);
// }

// Array<double> FFOps_Basic_neu::getGain(Individual& Ind){
//   Array<double> gain(getNNeurons(Ind));
//   for (unsigned i=0; i<gain.dim(0); i++)
//     gain(i) = ((ChromosomeT<double>&) (Ind[3]))[i];
//   return gain;
// }

Array<double> FFOps_Basic_neu::getGain(Individual& Ind){
  unsigned nIn  = getNIn(Ind);
  unsigned nHid = getNHid(Ind);
  unsigned nOut = getNOut(Ind);
  unsigned N    = nIn+nHid+nOut;
  Array<int> con  = getConnections(Ind,0);
  Array<double> wei  = getWeights(Ind,0);

  Array<double> gain_h(N);
  for (unsigned j=0; j<nIn; j++)
    gain_h(j)=0;

  for (unsigned j=nIn;j<N;j++){
    gain_h(j) = 0;
    for (unsigned i=0; i<j; i++)
      if (con(j,i)!=0)
        gain_h(j) += wei(j,i)*wei(j,i);
    
    gain_h(j) = sqrt(gain_h(j));
  }

  return gain_h;
}


void FFOps_Basic_neu::setActivation(Individual& Ind, Array<double> activation){
  if(activation.dim(0) != getNNeurons(Ind)){
    cerr << "error setActivation: dimension of activation vector (" << activation.dim(0) 
	 << ") and number of neurons (" << getNNeurons(Ind) << ") do not coincide." << endl;
    exit(EXIT_FAILURE);
  }

  ((ChromosomeT<double>&) (Ind[3])).resize(getNNeurons(Ind));
  for (unsigned i=0; i<activation.dim(0); i++)
    ((ChromosomeT<double>&) (Ind[3]))[i] = activation(i);
}

Array<double> FFOps_Basic_neu::getActivation(Individual& Ind){
  Array<double> activation(getNNeurons(Ind));
  for (unsigned i=0; i<activation.dim(0); i++)
    activation(i) = ((ChromosomeT<double>&) (Ind[3]))[i];
  return activation;
}

Array<int> FFOps_Basic_neu::getConnectionsWithoutBias(Individual& Ind){
  Array<int> cR = getConnectionsReClaM(Ind);
  Array<int> CWB(cR.dim(0),cR.dim(1)-1);
  for (unsigned i=0; i<CWB.dim(0); i++)
    for (unsigned j=0; j<CWB.dim(1); j++)
      CWB(i,j) = cR(i,j);
  return CWB;
}

Array<int> FFOps_Basic_neu::getConnectionsReClaM(Individual& Ind){
  Array<int> c = getConnections(Ind, 0);
  Array<int> cR(c.dim(0),c.dim(1)+1);
  cR = 0;
  for (unsigned z=0; z<c.dim(0); z++){
    for (unsigned s=0; s<z; s++)
      cR(z,s) = c(z,s);
    cR(z,cR.dim(1)-1) = c(z,z);
  }
  return cR;
}

Array<double> FFOps_Basic_neu::getWeightsReClaM(Individual& Ind){
  Array<double> w = getWeights(Ind, 0);
  Array<double> wR(w.dim(0),w.dim(1)+1);
  wR = 0;
  for (unsigned z=0; z<w.dim(0); z++){
    for (unsigned s=0; s<z; s++)
      wR(z,s) = w(z,s);
    wR(z,wR.dim(1)-1) = w(z,z);
  }
  return wR;
}

void FFOps_Basic_neu::setWeightsReClaM(Individual& Ind, Array<double> wR){
  if ( (wR.ndim()!=2) || (wR.dim(0)+1!=wR.dim(1)) ){
    cerr << "error in FFOps_Basic_neu::setWeightsReClaM: Dimensions do not fit " 
	 << "(" << wR.dim(0) << "," << wR.dim(1) << ")!" << endl;
    exit(EXIT_FAILURE);
  }else{
    Array<double> w(wR.dim(0),wR.dim(1)-1);
    w = 0;
    for (unsigned z=0; z<w.dim(0); z++){
      for (unsigned s=0; s<z; s++)
	w(z,s) = wR(z,s);
      w(z,z) = wR(z,wR.dim(1)-1);
    }
    setStructure(Ind, w, 0);
  }
}

void FFOps_Basic_neu::setConnectionsReClaM(Individual& Ind, Array<int> cR){
  if ((cR.dim(0)==cR.dim(1)) || (cR.dim(0)==cR.dim(1)-1)){
    Array<int> c(cR.dim(0),cR.dim(0));
    c = 0;
    for (unsigned z=0; z<c.dim(0); z++){
      for (unsigned s=0; s<z; s++)
	c(z,s) = cR(z,s);
      if (cR.dim(0)==cR.dim(1))
	if (z>=getNIn(Ind))
	  c(z,z)=1;
      else
	c(z,z) = cR(z,cR.dim(1)-1);
    }
    setStructure(Ind, c, 0);
  }
  else{
    cerr << "error in FFOps_Basic_neu::setConnectionsReClaM: dimensions do not fit" << endl; 
    exit(EXIT_FAILURE);
  }
}

void FFOps_Basic_neu::printNet(Individual& Ind, string filename){
  Array<double> wei = getWeightsReClaM(Ind);
  Array<int>    con = getConnectionsReClaM(Ind);
  unsigned      nIn = getNIn(Ind);
  unsigned      nOut = getNOut(Ind);
  FILE *fp=fopen(filename.c_str(),"wt");
  fprintf(fp,"%u %u\n",nIn, nOut);

  for (unsigned i=0; i<con.dim(0); i++){
    for (unsigned j=0; j<con.dim(1); j++){
      fprintf(fp,"%1d",con(i,j));
      if (j<con.dim(1)-1) fprintf(fp," ");
    }
    fprintf(fp,"\n");
  }

  fprintf(fp,"\n");

  for (unsigned i=0; i<wei.dim(0); i++){
    for (unsigned j=0; j<wei.dim(1); j++){
      fprintf(fp,"%10.8e",wei(i,j));
      if (j<wei.dim(1)-1) fprintf(fp," ");    
    }
    fprintf(fp,"\n");
  }

  fclose(fp);
}

void FFOps_Basic_neu::printNetWithoutBias(Individual& Ind, string filename){
  Array<double> wei = getWeightsReClaM(Ind);
  Array<int>    con = getConnectionsWithoutBias(Ind);
  unsigned      nIn = getNIn(Ind);
  unsigned      nOut = getNOut(Ind);
  FILE *fp=fopen(filename.c_str(),"wt");
  fprintf(fp,"%u %u\n",nIn, nOut);

  for (unsigned i=0; i<con.dim(0); i++){
    for (unsigned j=0; j<con.dim(1); j++){
      fprintf(fp,"%1d",con(i,j));
      if (j<con.dim(1)-1) fprintf(fp," ");
    }
    fprintf(fp,"\n");
  }

  fprintf(fp,"\n");

  for (unsigned i=0; i<wei.dim(0); i++){
    for (unsigned j=0; j<wei.dim(1); j++){
      fprintf(fp,"%10.8e",wei(i,j));
      if (j<wei.dim(1)-1) fprintf(fp," ");    
    }
    fprintf(fp,"\n");
  }

  fclose(fp);
}

//#pragma implementation
#include "./Ops_Basic.h"

//////////////////////////////////////////////////////////////
//
//  Initialize
//
//////////////////////////////////////////////////////////////
Individual Ops_Basic::createEmpty(const unsigned nIn, const unsigned nHid, const unsigned nOut, unsigned nDelay){
  unsigned N = nIn + nHid + nOut;
  unsigned laenge = N*N;
  Individual Ind(ChromosomeT <unsigned>    (4));
  for (unsigned i=0; i<nDelay; i++){
    ChromosomeT<int> con(laenge);
    for (unsigned j=0; j<con.size(); j++)
      con[j]=0;
    Ind.append(con);
    ChromosomeT<double> wei(laenge);
    for (unsigned j=0; j<con.size(); j++)
      wei[j]=0.0;
    Ind.append(wei);
  }

  ((ChromosomeT<unsigned>&) (Ind[0]))[0] = nIn;
  ((ChromosomeT<unsigned>&) (Ind[0]))[1] = nHid;
  ((ChromosomeT<unsigned>&) (Ind[0]))[2] = nOut;
  ((ChromosomeT<unsigned>&) (Ind[0]))[3] = nDelay;

  return Ind;
}

void Ops_Basic::initWeights(Individual& Ind, double low, double high){
  Array<int>    c = getConnections(Ind);
  Array<double> w = getWeights(Ind);
  
  for (unsigned i=0; i<c.dim(0); i++)
    for (unsigned j=0; j< c.dim(1); j++)
      for (unsigned k=0; k<c.dim(2); k++)
	if (c(i,j,k)!=0)
	  w(i,j,k) = Rng::uni(low,high);
  
  setStructure(Ind, w);
}


//////////////////////////////////////////////////////////////
//
// Kleinzeug
//
//////////////////////////////////////////////////////////////

unsigned Ops_Basic::getNIn(Individual& Ind){
  return ((ChromosomeT<unsigned>&) (Ind[0]))[0];}

unsigned Ops_Basic::getNHid(Individual& Ind){
  return ((ChromosomeT<unsigned>&) (Ind[0]))[1];}

unsigned Ops_Basic::getNOut(Individual& Ind){
  return ((ChromosomeT<unsigned>&) (Ind[0]))[2];}

unsigned Ops_Basic::getNDelay(Individual& Ind){
  return ((ChromosomeT<unsigned>&) (Ind[0]))[3];}

unsigned Ops_Basic::getNNeurons(Individual& Ind){
  return getNIn(Ind)+getNHid(Ind)+getNOut(Ind);}

void Ops_Basic::setNHid(Individual& Ind, unsigned nHid){
  ((ChromosomeT<unsigned>&) (Ind[0]))[1] = nHid;
  
  unsigned N = getNIn(Ind) + nHid + getNOut(Ind);
  unsigned laenge = N*N;
  unsigned nDelay = getNDelay(Ind);
  for (unsigned i=1; i<2*nDelay+1; i+=2){
    ((ChromosomeT<int>&) (Ind[i])).resize(laenge);
    ((ChromosomeT<double>&) (Ind[i+1])).resize(laenge);
  }
}

//////////////////////////////////////////////////////////////
//
// getConnections
//
//////////////////////////////////////////////////////////////

Array<int> Ops_Basic::getConnections(Individual& Ind){
  unsigned nDelay = getNDelay(Ind);
  Array<int> c(nDelay,getNNeurons(Ind),getNNeurons(Ind));  
  for (unsigned i=0; i<nDelay; i++)
    c[i] = getConnections(Ind, i);
  
  return c;
}

Array<int> Ops_Basic::getConnections(Individual& Ind, unsigned entry){
  unsigned N=getNNeurons(Ind);
  Array<int> connections(N,N);
  entry = 2*entry+1;
  
  unsigned count=0;
  for (unsigned i=0; i<N; i++)
    for (unsigned j=0; j<N; j++){
      connections(i,j) = ((ChromosomeT<int>&) (Ind[entry]))[count];
      count++;
    }

  return connections;
}

unsigned Ops_Basic::getNParameter(Individual& Ind){
  Array<int> c;
  c = getConnections(Ind);
  unsigned nParams = 0;
  for (unsigned i=0; i<c.dim(0); i++)
    for (unsigned j=0; j<c.dim(1); j++)
      for (unsigned k=0; k<c.dim(2); k++)
	if (c(i,j,k) != 0)
	  nParams++;
  return nParams;
} 

//////////////////////////////////////////////////////////////
//
// getWeights
//
//////////////////////////////////////////////////////////////

Array<double> Ops_Basic::getWeights(Individual& Ind, unsigned entry){
  unsigned N=getNNeurons(Ind);
  Array<double> weights(N,N);
  entry = 2*(entry+1);
  
  unsigned count=0;
  for (unsigned i=0; i<N; i++)
    for (unsigned j=0; j<N; j++){
      weights(i,j) = ((ChromosomeT<double>&) (Ind[entry]))[count];
      count++;
    }

  return weights;
}

Array<double> Ops_Basic::getWeights(Individual& Ind){
  unsigned nDelay = getNDelay(Ind);
  Array<double> w(nDelay,getNNeurons(Ind),getNNeurons(Ind));
  for (unsigned i=0; i<nDelay; i++)
    w[i] = getWeights(Ind, i);
  
  return w;
}

void Ops_Basic::printNet(Individual& Ind, string filename){
  unsigned nIn = getNIn(Ind);
  unsigned nOut = getNOut(Ind);
  FILE *fp=fopen(filename.c_str(),"wt");
  fprintf(fp,"%u %u\n",nIn, nOut);

  unsigned nDelay = getNDelay(Ind);
  for (unsigned k=0; k<nDelay; k++){  
    Array<int> con = getConnections(Ind, k);
    for (unsigned i=0; i<con.dim(0); i++){
      for (unsigned j=0; j<con.dim(1); j++){
	fprintf(fp,"%1d",con(i,j));
	if (j<con.dim(1)-1) fprintf(fp," ");
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }

  for (unsigned k=0; k<nDelay; k++){  
    Array<double> wei = getWeights(Ind, k);
    for (unsigned i=0; i<wei.dim(0); i++){
      for (unsigned j=0; j<wei.dim(1); j++){
	fprintf(fp,"%10.8e",wei(i,j));
	if (j<wei.dim(1)-1) fprintf(fp," ");
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

//////////////////////////////////////////////////////////////
//
// setStructure
//
//////////////////////////////////////////////////////////////

void Ops_Basic::setStructure(Individual& Ind, const Array<double>& w, unsigned entry){
  if ( (w.ndim()!=2) || (w.dim(0)!=w.dim(1)) ){
    cerr << "error in Ops_Basic::setStructure: Dimensions do not fit." << endl;
    exit(EXIT_FAILURE);
  }else{
    unsigned N=w.dim(0);
    entry = 2*(entry+1);
    
    unsigned count=0;
    for (unsigned i=0; i<N; i++)
      for (unsigned j=0; j<N; j++){
	((ChromosomeT<double>&) (Ind[entry]))[count]=w(i,j);
	count++;
      }
  }
}

void Ops_Basic::setStructure(Individual& Ind, const Array<int>& c, unsigned entry){
  if ( (c.ndim()!=2) || (c.dim(0)!=c.dim(1)) ){
    cerr << "error in Ops_Basic::setStructure: Dimensions do not fit." << endl;
    exit(EXIT_FAILURE);
  }else{
    unsigned N=c.dim(0);
    entry = 2*entry+1;
    
    unsigned count=0;
    for (unsigned i=0; i<N; i++)
      for (unsigned j=0; j<N; j++){
	((ChromosomeT<int>&) (Ind[entry]))[count]=c(i,j);
	count++;
      }
  }
}

void Ops_Basic::setStructure(Individual& Ind, const Array<double>& w){
  if ((w.ndim()!=3) || (w.dim(0)!=getNDelay(Ind)) || (w.dim(1)!=w.dim(2))){
    cerr << "error in Ops_Basic::setStructure: Dimensions do not fit." << endl;
    exit(EXIT_FAILURE);
  }else{
    unsigned nDelay = getNDelay(Ind);
    for (unsigned i=0; i<nDelay; i++)
      setStructure(Ind,w[i],i);
  }
}

void Ops_Basic::setStructure(Individual& Ind, const Array<int>& c){
  if ((c.ndim()!=3) || (c.dim(0)!=getNDelay(Ind)) || (c.dim(1)!=c.dim(2))){
    cerr << "error in Ops_Basic::setStructure: Dimensions do not fit." << endl;
    exit(EXIT_FAILURE);
  }else{
    unsigned nDelay = getNDelay(Ind);
    for (unsigned i=0; i<nDelay; i++)
      setStructure(Ind,c[i],i);
  }
}

void Ops_Basic::delNeuron(Individual& Ind, unsigned N){
  unsigned nIn = getNIn(Ind);
  if (N < getNHid(Ind)){
    unsigned nDelay = getNDelay(Ind);
    N += nIn;
    Array<int>    c = getConnections(Ind),
      C(nDelay,c.dim(1)-1,c.dim(2)-1);
    Array<double> w = getWeights(Ind),
      W(nDelay,w.dim(1)-1,w.dim(2)-1);

    C = 0;
    W = 0;

    for (unsigned i=0; i<nDelay; i++){
      for (unsigned z=0; z<N; z++){
	for (unsigned s=0; s<N; s++){
	  C(i,z,s) = c(i,z,s);
	  W(i,z,s) = w(i,z,s);
	}
	for (unsigned s=N; s<C.dim(2); s++){
	  C(i,z,s) = c(i,z,s+1);
	  W(i,z,s) = w(i,z,s+1);
	}
      }
    
      for (unsigned z=N; z<C.dim(1); z++){
	for (unsigned s=0; s<N; s++){
	  C(i,z,s) = c(i,z+1,s);
	  W(i,z,s) = w(i,z+1,s);
	}
	for (unsigned s=N; s<C.dim(2); s++){
	  C(i,z,s) = c(i,z+1,s+1);
	  W(i,z,s) = w(i,z+1,s+1);
	}
      }
    }
    
    //    Ind = createEmpty(nIn,getNHid(Ind)-1,getNOut(Ind),nDelay);
    setNHid(Ind, getNHid(Ind)-1);
    setStructure(Ind,C);
    setStructure(Ind,W);
  }
}

void Ops_Basic::insNeuron(Individual& Ind, unsigned N){
  if (N <= getNHid(Ind)){
    N += getNIn(Ind);
    Array<int>    c = getConnections(Ind),
      C(c.dim(0),c.dim(1)+1,c.dim(2)+1);
    Array<double>    w = getWeights(Ind),
      W(w.dim(0),w.dim(1)+1,w.dim(2)+1);

    C=0;
    W=0;

    for (unsigned i=0; i<getNDelay(Ind); i++){
      for (unsigned z=0; z<N; z++){
	for (unsigned s=0; s<N; s++){
	  C(i,z,s) = c(i,z,s);
	  W(i,z,s) = w(i,z,s);
	}
	for (unsigned s=N; s<c.dim(2); s++){
	  C(i,z,s+1) = c(i,z,s);
	  W(i,z,s+1) = w(i,z,s);
	}
      }
    
      for (unsigned z=N; z<c.dim(1); z++){
	for (unsigned s=0; s<N; s++){
	  C(i,z+1,s) = c(i,z,s);
	  W(i,z+1,s) = w(i,z,s);
	}
	for (unsigned s=N; s<c.dim(2); s++){
	  C(i,z+1,s+1) = c(i,z,s);
	  W(i,z+1,s+1) = w(i,z,s);
	}
      }
    }
    
    //    Ind = createEmpty(getNIn(Ind),getNHid(Ind)+1,getNOut(Ind),getNDelay(Ind));
    setNHid(Ind, getNHid(Ind)+1);
    setStructure(Ind,C);
    setStructure(Ind,W);
  }
}

//////////////////////////////////////////////////////////////
//
// jog weights
//
//////////////////////////////////////////////////////////////
void Ops_Basic::jogWeights(Individual& Ind, const double stdv, double prob){
  Array<int> c = getConnections(Ind);
  Array<double> w = getWeights(Ind);
  double stdv2 = stdv*stdv;
  unsigned neurons = c.dim(1);
  for (unsigned i=0; i<c.dim(0); i++)
    for (unsigned j=0; j<neurons; j++)    
      for (unsigned k=0; k<neurons; k++)
	if ((c(i,j,k) != 0) && (Rng::uni(0,1)<prob))
	  w(i,j,k) += Rng::gauss(0,stdv2);

  setStructure(Ind, w);
}

//////////////////////////////////////////////////////////////
//
// Gluecksrad
//
//////////////////////////////////////////////////////////////
unsigned Ops_Basic::wheelOfFortune(Array<double>& prob){
  for (unsigned i=1; i<prob.dim(0); i++)
    prob(i) += prob(i-1);
  
  double number = Rng::uni(0,prob(prob.dim(0)-1));
  unsigned pos = 0;
  while (number > prob(pos))
    pos ++;
  
  return pos;
}

//////////////////////////////////////////////////////////////
//
// pvm-routines
//
//////////////////////////////////////////////////////////////

/*
int Ops_Basic::pvm_pkIndividual(Individual& Ind){
  unsigned nChrom = Ind.size();
  
  // sende Laenge des Individuums
  int *d = new int[1];
  d[0]=nChrom;
  pvm_pkint(d,1,1);
  delete[] d;

  // sende einzelne Chromosomen
  for (unsigned i=0; i<nChrom; i++){
    // sende Infos
    ChromosomeT<double> doubleType(1);
    ChromosomeT<int> intType(1); 
    ChromosomeT<unsigned> unsignedType(1); 
    int typ=0;
    if ( Ind[i].sameType(doubleType) )
      typ = 0;
    if ( Ind[i].sameType(intType) )
      typ = 1;
    if ( Ind[i].sameType(unsignedType) )
      typ = 2;
    unsigned laenge = Ind[i].size();
    int *d = new int[2];
    d[0]=typ;
    d[1]=(int)laenge;
    pvm_pkint(d,2,1);
    delete[] d;   

    // sende Werte
    if ( typ==0 ){ 
      vector<double> v=(ChromosomeT<double>&)(Ind[i]);
      double *werte = new double[laenge];
      for (unsigned j=0; j<laenge; j++)
	werte[j]=((ChromosomeT<double>&)(Ind[i]))[j];
      pvm_pkdouble(werte,(int)laenge,1);
      delete[] werte;
    }
 
    if ( typ==1 ){ 
      int *werte = new int[laenge];
      for (unsigned j=0; j<laenge; j++)
	werte[j]=((ChromosomeT<int>&)(Ind[i]))[j];
      pvm_pkint(werte,(int)laenge,1);
      delete[] werte;
    }
 
    if ( typ==2 ){ 
      int *werte = new int[laenge];
      for (unsigned j=0; j<laenge; j++)
	werte[j]=((ChromosomeT<unsigned>&)(Ind[i]))[j];
      pvm_pkint(werte,(int)laenge,1);
      delete[] werte;
    }
  }

  // Fitness des Individuums;
  double *f = new double[1];
  f[0] = Ind.fitnessValue();
  pvm_pkdouble(f,1,1);
  delete[] f;

  //  cout << typeInfo(Ind(0)) << endl;
  //   vector<unsigned> v=(ChromosomeT<unsigned>&)(Ind[0]);
//   unsigned size = v.size();
//   int *d = new int[size];
//   for (unsigned i=0; i<size; i++)
//     d[i] = (int)v[i];
//   pvm_pkint(d,size,1);
//   delete[] d;

//   double *f = new double[1];
//   f[0] = Ind.fitnessValue();
//   pvm_pkdouble(f,1,1);
//   delete[] f;

//   for (unsigned i=0; i<getNDelay(Ind); i++){
//     vector<int> con = (ChromosomeT<int>&) (Ind[2*i+1]);
//     size = con.size();
//     int *c = new int[size];
//     for (unsigned j=0; j<size; j++)
//       c[j] = con[j];
//     pvm_pkint(c,size,1);
//     delete[] c;
    
//     vector<double> wei = (ChromosomeT<double>&) (Ind[2*i+2]);
//     size = wei.size();
//     double *w = new double[size];
//     for (unsigned j=0; j<size; j++)
//       w[j] = wei[j];
//     pvm_pkdouble(w,size,1);
//     delete[] w;
//   }
  
  return 1;
}

int Ops_Basic::pvm_upkIndividual(Individual& Ind){
  // Laenge des Individuums
  int *d = new int[1];
  pvm_upkint(d,1,1);
  unsigned nChrom = d[0];
  delete[] d;

  Individual Ind_help;
  for (unsigned i=0; i<nChrom; i++){
    // Info ueber Chromosome
    int *info = new int[2];
    pvm_upkint(info,2,1);
    int typ = info[0];
    unsigned laenge = (unsigned) info[1];
    delete[] info;

    // lies eigentliches Chromosome;
    if (typ==0){
      double *werte = new double[laenge];
      pvm_upkdouble(werte,laenge,1);
      ChromosomeT<double> Chrom(laenge);
      for (unsigned j=0; j<laenge; j++)
	Chrom[j] = werte[j];
      delete [] werte;
      Ind_help.append(Chrom);     
    }
    if (typ==1){
      int *werte = new int[laenge];
      pvm_upkint(werte,laenge,1);
      ChromosomeT<int> Chrom(laenge);
      for (unsigned j=0; j<laenge; j++)
	Chrom[j] = werte[j];
      delete [] werte;
      Ind_help.append(Chrom);     
    }
    if (typ==2){
      int *werte = new int[laenge];
      pvm_upkint(werte,laenge,1);
      ChromosomeT<unsigned> Chrom(laenge);
      for (unsigned j=0; j<laenge; j++)
	Chrom[j] = (unsigned) werte[j];
      delete [] werte;
      Ind_help.append(Chrom);     
    }
  }

  // Fitness des Individuums
  double *f = new double[1];
  pvm_upkdouble(f,1,1);
  Ind_help.setFitness(f[0]);
  delete[] f;

//   unsigned size = 4;
//   int *d = new int[size];
//   pvm_upkint(d,size,1);
//   Individual Ind_help;
//   ChromosomeT<unsigned> Chrom1(size);
//   for (unsigned i=0; i<size; i++)
//     Chrom1[i] = (unsigned) d[i];
//   Ind_help.append(Chrom1);
//   size = (d[0]+d[1]+d[2])*(d[0]+d[1]+d[2]);
//   unsigned nDelay = d[3];
//   delete[] d;

//   double *f = new double[1];
//   pvm_upkdouble(f,1,1);
//   Ind_help.setFitness(f[0]);
//   delete[] f;
  
//   for (unsigned i=0; i<nDelay; i++){
//     int *c = new int[size];
//     pvm_upkint(c,size,1);
//     ChromosomeT<int> ChromC(size);
//     for (unsigned j=0; j<size; j++)
//       ChromC[j] = c[j];
//     delete[] c;
//     Ind_help.append(ChromC);

//     double *w = new double[size];
//     pvm_upkdouble(w,size,1);
//     ChromosomeT<double> ChromW(size);
//     for (unsigned j=0; j<size; j++)
//       ChromW[j] = w[j];
//     delete[] w;
//     Ind_help.append(ChromW);
//   }

  Ind=Ind_help;
  return 1;
}
*/




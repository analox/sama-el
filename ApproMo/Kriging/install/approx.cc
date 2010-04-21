#include "approx.h"

approx::approx() {
  p = 0;
  n = 0;
  X.newsize(1,1);
  X = 0;
  y.newsize(1);
  y = 0;
} // end constructor

approx::~approx() {}

void approx::notseterror(){
  cout << "ERROR: Value not set." << endl;
  assert(0==1);  // FAILS!!!
}

long approx::getp() {
   if ( !ValuesSet )
        notseterror();
   return p;
}

long approx::getn() {
   if ( !ValuesSet )
        notseterror();
   return n;
}

Matrix<double> approx::getPoints(){
    if ( !ValuesSet )
        notseterror();
    return (X);
}

Vector<double> approx::getValues(){
    if ( !ValuesSet )
        notseterror();

    return y;
}


double reportsample( approx& Func, Vector<double>& x, ostream& outf) {
  long p = x.size();
  double result= 0;
  long i;

  for( i = 0; i < p; i++ ) {
    outf << x[i] << ' ';
  } // end for

  result = Func.evalf(x);
  outf << result << endl;
  return result;
} // end reportsample()
   
void scanfunction( approx& Func, ostream& outf, long pointsperaxis,
                   Vector<double>& lower, Vector<double>& upper )
  // Last Modified 7/14/99
{ 
  long p = Func.getp();
  Vector<double> x(p);
  Vector<double> increment(p);
  long j;
  double temp = 1.0/pointsperaxis;

  if( p == 2) {  // use the faster version if possible
    scanplane( Func, outf, pointsperaxis, lower, upper);
    return;
  } // end if


  x = lower;
  increment = upper-lower;
  increment = increment * temp;

  while( x[0] < (upper[0]+increment[0]/2.0) ) {
    reportsample( Func, x, outf ); 
    x[p-1] += increment[p-1];
    for( j = p-1; (j > 0)&&(x[j] > (upper[j]+increment[j]/2.0)); j--) {
      x[j] = lower[j];
      x[j-1] += increment[j];
      if( j-1 == 0)
	outf << endl;
    } // end for
  } // end while
   
  return;
} // end scanfunction

void scanfunction( approx& Func, ostream& outf, Vector<double>& lower,
                   Vector<double>& upper ) {
  scanfunction(Func, outf, DEFAULTPOINTS, lower, upper);
}
// A default version which uses DEFAULTPOINTS pointsperaxis.  You can modify this default
// by changing the value of DEFAULTPOINTS in the approx.h file.


 void scanplane( approx& Func, ostream& outf, long pointsperaxis,
                 Vector<double>& lower, Vector<double>& upper )  // Last Modified 7/15/99
{ 
  double value = 0, stepe, stepn;
  Vector<double> x(2);

  stepe = (upper[0]-lower[0])/pointsperaxis;
  stepn = (upper[1]-lower[1])/pointsperaxis;

  for( x[0]=lower[0]; x[0]<(upper[0]+stepe/2.0); x[0]+= stepe ) {
    for( x[1]=lower[1]; x[1]<(upper[1]+stepn/2.0); x[1]+=stepn ) {
      outf << x[0] << ' ' << x[1] << ' ';
      value = Func.evalf(x);
      outf << value << endl;
    } // end for
    outf << endl;
  } // end for
  return;
} // end scanplane


void analyzevalues( approx& Func, ostream& outf, Vector<double>& lower, Vector<double>& upper ) {
  analyzevalues( Func, outf, 25, lower, upper);
  return;
}
/****
     Works much like scanfunction, except that it reports simple statistics on
     the points sampled to standard I/O, such as maximum value, minimum value,
     average value, and the approximate location of the 25th and 75th
     percentiles.

     Warning:  This function stores all sampled points in memory.  In higher
     dimensions, this may fill up the available memory.  Use at your own risk.
****/


void analyzevalues( approx& Func, ostream& outf, long pointsperaxis, Vector<double>& lower, Vector<double>& upper )
  /*Accepts a krigapprox object and two filenames.  Will create output files using
    the filenames.  filename1 will be the result of scanfunction().  filename 2
    will contain an analysis of these results.  While not exact, this analysis
    should provide the user with an idea of the range of values the function
    assumes.

  */
{
  double max, min, avg=0;
  long psize = Func.getp(), i, numpoints = 0, numpoints2 = 0, count = 0;
  Vector<double> current(1), maxvec(1), minvec(1);

  long p =psize;
  Vector<double> x(p);
  Vector<double> increment(p);
  long j;
  double temp = 1.0/((double)pointsperaxis);
  Matrix<double> points (1,1);

  x = lower;
  increment = upper-lower;
  increment = increment * temp;
  temp = pow((double)(pointsperaxis+2.0), (double)p);
  points.newsize((long)temp,p+1);
  
  while(( x[0] < (upper[0]+increment[0]/2.0) )&& (count<(long)temp)) {

    for(long k = 0; k <p; k++) {
      points[count][k] = x[k];
      outf << x[k] << ' ';
    } // end for
    points[count][p] = Func.evalf(x);
    outf << points[count][p] << endl;
    count++;
    
    x[p-1] += increment[p-1];
    for( j = p-1; (j > 0)&&(x[j] > (upper[j]+increment[j]/2.0)); j--) {
      x[j] = lower[j];
      x[j-1] += increment[j];
      if( j-1 == 0)
	outf << endl;
    } // end for
  } // end while
  temp = count;

  // Process results of scanfunction()
  current.newsize(psize+1);
  maxvec.newsize(psize+1);
  minvec.newsize(psize+1);  

  count = 1;
  current = points.row(count);
  count++;
  avg = avg + current[psize];
  numpoints++;
  maxvec = current;
  minvec = current;
  max = maxvec[psize];
  min = minvec[psize];
  
  while( count < temp ) {
    current = points.row(count);
    count++;
    
    avg = avg+current[psize];
    numpoints++;
    if( current[psize] > max ) {
      maxvec = current;
      max = maxvec[psize];
    } else if ( current[psize] < min ) {
       minvec = current;
       min = minvec[psize];
    } // end elseif
  } // end while

  avg = avg/numpoints;

  cout <<  "\nMaximum value of " << max << " found at (" << maxvec[0];
  for( i=1; i<psize; i++)
    cout << ", " << maxvec[i];
  cout << ")" << endl;
  cout <<  "Minimum value of " << min << " found at (" << minvec[0];
  for( i=1; i<psize; i++)
    cout << ", " << minvec[i];
  cout << ")" << endl;
  
  cout << "The average value of the points sampled was " << avg << endl;
  
  count = 1;

  numpoints = 0;
  numpoints2 = 0;
  max = 0;
  min = 0;
  
  current = points.row(count);
  count++;
  if( current[psize] >= avg) {
    max += current[psize];
    numpoints++;
  } else {
    min += current[psize];
    numpoints2++;
  }
    
  while( count < temp) {
    current = points.row(count);
    count++;
    if( current[psize] >= avg) {
      max += current[psize];
      numpoints++;
    } else {
      min += current[psize];
      numpoints2++;
    } // end else
  } // end while

  max = max/numpoints;
  min = min/numpoints2;

  cout << "\n25th percentile is approximately " << min << endl;
  cout << "75th percentile is approximately " << max << endl;

  return;
}  // end analyzevalues()



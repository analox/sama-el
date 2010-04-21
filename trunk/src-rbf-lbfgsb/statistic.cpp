#include<iostream.h>
#include<fstream.h>
#include<stdlib.h>

int compare_numb(const void *num1, const void *num2)
{
  double a, b;
	a = (double)(*(double*)num1);
 	b = (double)(*(double*)num2);





  if (a < b) return -1;
  if (a == b) return  0;
  if (a > b) return  1;
}




int main(int argc, char** argv) {

	unsigned i, j, nFile, nEval;
	
	if(argc<3) {
		cout << "Please put arg1 = nFile, arg2 = nEval" << endl;
		exit(1);
	}

	nFile = atoi(argv[1]);
	nEval = atoi(argv[2]);
	double index, value[nFile][nEval], averaged[nEval];
	for(i=0; i<nFile; i++) {
		char argu[30];
		sprintf(argu,"logs/bestFile-%d.dat",i+1);
		
		ifstream input(argu);
		for(j=0;j<nEval;j++) {
			input >> index >> value[i][j];
		
		}
		input.close();
	}
	
	ofstream output("statistic.dat", ios::trunc);
	output << "#index\t\tmean\t\tstdev\t\tmedian\t\tbest\t\tworst" << endl;
	for(i=0; i<nEval; i++) {
		double sum = 0;
		for(j=0; j<nFile; j++) {
			sum += value[j][i];
		}
		averaged[i] = sum/nFile;
		//output << i << "\t" << averaged[i] << endl;


		double a[nFile];
		for(j=0; j<nFile; j++) {
			a[j] = value[ j ][ i ];
		}

		qsort(
		        a, 			/* Pointer to elements		*/
        		nFile, 			/* Number of elements		*/
       	 		sizeof(double),  	/* size of one element.		*/
        		compare_numb		/* Pointer to comparison function */
       		);

		// output average, median, best, worst

		double stdev=0;
		for(j=0; j<nFile; j++) {
			stdev += (value[j][i]-averaged[i])*(value[j][i]-averaged[i]);
		}
		stdev /= (double)(nFile-1);
		stdev = sqrt( stdev );

		if(nFile%2 == 0 )output << i << "\t\t" << averaged[ i ] << "\t\t" << stdev << "\t\t" << (a[(nFile/2)-1] + a[nFile/2]) / 2.0 << "\t\t" << a[0] << "\t\t" << a[nFile-1] << endl;
		else output << i << "\t\t" << averaged[ i ] << "\t\t" << stdev << "\t\t" << a[ (int)(floor( nFile/2 )) ] << "\t\t" << a[0] << "\t\t" << a[nFile-1] << endl;
		
	}
	output.close();
	
	
	
	return 0;

}

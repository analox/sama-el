#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<math.h>
using namespace std;
int compare_numb(const void *num1, const void *num2)
{
  double a, b;
	a = (double)(*(double*)num1);
 	b = (double)(*(double*)num2);

  if (a < b) return -1;
  if (a == b) return  0;
  if (a > b) return  1;
}


double _index, value[50][300000], averaged[300000], maxEval[50];

int main(int argc, char** argv) {

	unsigned i, j, nFile, nEval;
	
	if(argc<3) {
		cout << "Please put arg1 = nFile, arg2 = nEval" << endl;
		exit(1);
	}

	nFile = atoi(argv[1]);
	nEval = atoi(argv[2]);
	
	cout << "nFile: " << nFile << ", nEval: " << nEval << endl;


	for(i=0; i<nFile; i++) {

		maxEval[i] = nEval;

		char argu[30];
		sprintf(argu,"logs/bestFile-%d.dat",i+1);
		
		ifstream input(argu);
		for(j=0;j<nEval;j++) {
			input >> _index >> value[i][j];
			if (input.eof()) {
				cout << "Stop at: " << j << endl;
				maxEval[i] = j;
				break;
			}
		
		}
		input.close();
	}
	
	ofstream output("statistic.dat", ios::trunc);
	output << "#index\t\tmean\t\tstdev\t\tmedian\t\tbest\t\tworst" << endl;

	double a[50];	
	for(i=0; i<nEval; i++) {
		double sum = 0;
		int fileCount = 0;
		for(j=0; j<nFile; j++) {
			if (i < maxEval[j]) {
				// Calculate sum of valid value
				sum += value[ j ][ i ];
				// Copy valid value to a
				a[fileCount] = value[ j ][ i ];
				fileCount++;
			}
		}

		if (fileCount == 0) 
		{
			cout << "Exceeding evaluation..." << endl;
			break;
		}
		else if (fileCount < nFile)
			cout << "File count: " << fileCount << endl;

		averaged[i] = sum/fileCount;
		//output << i << "\t" << averaged[i] << endl;

		qsort(
		        a, 			/* Pointer to elements		*/
        		fileCount,		/* Number of elements		*/
       	 		sizeof(double),  	/* size of one element.		*/
        		compare_numb		/* Pointer to comparison function */
       		);

		// output average, median, best, worst
		double stdev=0;
		for(j=0; j<fileCount; j++) {
			stdev += (a[j]-averaged[i])*(a[j]-averaged[i]);
		}
		stdev /= (double)(fileCount-1);
		stdev = sqrt( stdev );

		if(fileCount%2 == 0 )
			output << i << "\t\t" << averaged[ i ] << "\t\t" << stdev << "\t\t" << (a[(fileCount/2)-1] + a[fileCount/2]) / 2.0 << "\t\t" << a[0] << "\t\t" << a[fileCount-1] << endl;
		else 
			output << i << "\t\t" << averaged[ i ] << "\t\t" << stdev << "\t\t" << a[ (int)(floor( fileCount/2 )) ] << "\t\t" << a[0] << "\t\t" << a[fileCount-1] << endl;
		
	}
	output.close();
	
	return 0;

}

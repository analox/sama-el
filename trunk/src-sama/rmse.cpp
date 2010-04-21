#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;
#define MAX_GEN 50
#define MAX_MODEL 10

int main(int argc, char ** argv)
{
	char s[1000001];
  	FILE* f = fopen(argv[1], "r");
	int nGen = -1;
	int count = 0;
	int nModel = 3;
	double rmse[MAX_MODEL] = {0};
	int nS[MAX_MODEL] = {0};
	double fitStart, fitOptLoc, fitStartApprox, fitOptLocApprox;
	int modelFlag = 0;
	
	double res[MAX_GEN][MAX_MODEL] = {0};

	while(fgets(s, 1000000, f) != NULL)
	{
		string str(s);
		if (str.find("generation") != string::npos)
		{
			// restart counter
			if (str.find("generation 1 =") != string::npos) 
			{
				/* DEBUG 
				for (int i = 0; i < nGen; i++)
				{
					cout << (i+1) << ", ";
					for (int j = 0; j < nModel; j++) 
					{
						if (j < (nModel-1)) cout << res[i][j] << ", ";
						else cout << res[i][j] << endl;
					}
				}
				
				cout << "-----------------" << endl << str ;
				*/

				// To prevent interuption when building model in the end 
				// of last run - 1st RUN SURE NO LEARNING is PERFORMED
				for (int i = 0; i < nModel; i++) {
					rmse[i] = 0;
					nS[i] = 0;
				}

				if (nGen < count) nGen = count;
				count = 0;
			}
			//cout << str;
			//cout << "Gen " << count << endl;			
			double tmp = 0;
			for (int i = 0; i < nModel; i++) 
			{
				// cout << "E: " << rmse[i] << ", nS: " << nS[i] << endl;
				// calculate actual RMSE
				if (nS[i] > 0) rmse[i] = sqrt(rmse[i]/nS[i]);
				tmp += rmse[i];
			}
			for (int i = 0; i < nModel; i++) 
			{		
				// record accumulate normalized RMSE 
				if (tmp > 0)
					res[count][i] += rmse[i]/tmp;
				else 
					res[count][i] += rmse[i];
			}
			//cout << rmse[0] << ", " << rmse[1] << ", " << rmse[2] << endl;
			for (int i = 0; i < nModel; i++) {
				rmse[i] = 0;
				nS[i] = 0;
			}
			count++;
		}
		else if (str.find("BUILDING RBF") != string::npos)
		{
			//cout << str;
			modelFlag = 0;
		}
		else if (str.find("BUILDING PR") != string::npos)
		{
			//cout << str;
			modelFlag = 1;
		}
		else if (str.find("BUILDING KRIG") != string::npos)
		{
			//cout << str;
			modelFlag = 2;
		}
		else if (str.find("fitStart") != string::npos)
		{
			double tmp;
			sscanf(s, "fitStart = %lf", &tmp);
			fitStart = tmp;
			//cout << "fitStart = " << fitStart << endl;

			fgets(s, 1000000, f);
			sscanf(s, "fitOptLoc = %lf", &tmp);
			fitOptLoc = tmp;
			//cout << "fitOptLoc = " << fitOptLoc << endl;

			fgets(s, 1000000, f);
			sscanf(s, "fitStartApprox = %lf", &tmp);
			fitStartApprox = tmp;
			//cout << "fitStartApprox = " << fitStartApprox << endl;

			fgets(s, 1000000, f);
			sscanf(s, "fitOptLocApprox = %lf", &tmp);
			fitOptLocApprox = tmp;
			//cout << "fitOptLocApprox = " << fitOptLocApprox << endl;

			// Calculate for modelFlag
			double e = 0;
			e = fabs(fitStart - fitStartApprox);
			//cout << "Error: "<< e << endl;
			rmse[modelFlag] += e*e;
			nS[modelFlag] += 1;

			e = fabs(fitOptLoc - fitOptLocApprox);
			rmse[modelFlag] += e*e;
			nS[modelFlag] += 1;
			
		}
	}
	fclose(f);

	if (nGen < count) nGen = count;
	// Normalize RMSE
	for (int i = 0; i < nGen; i++)
	{
		//cout << "Gen " << i << ": ";
		double tmp = 0;
		for (int j = 0; j < nModel; j++) tmp += res[i][j];
		for (int j = 0; j < nModel; j++) 
		{	
			if (tmp > 0) res[i][j] = res[i][j]/tmp;
		}
		for (int j = 0; j < nModel; j++) 
		{
			if (j < (nModel-1)) cout << res[i][j] << ", ";
			else cout << res[i][j] << endl;
		}
	}
}

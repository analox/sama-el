#include "oss2.cpp"
#include "../tools/utilities.h"
int main()
{
	cout << "==========================================\n";
	cout << "Loading bounds..." << endl;
	ifstream f("usedOSS2param.param");

	int nparam=0;
	int i, j, k;
	double params[50];
	while(f >> params[nparam])
	{
		//cout << "  " << params[nparam] << endl;
		nparam++;
	}
	cout << "Reading inputs..." << endl;
	read_input();
	cout << "Calculate energy..." << endl;
	FEnergy fE;
	long long startT = getTime();
	for (j = 0; j < 3; j++) {
		for (i = 0; i < 100; i++) {
			cout << i << " " << data[i].nAtoms << " " << data[i].nO << " ";
			double energy;
			double grad[100][3];
			fE.energy_grad(data[i].r, data[i].nAtoms, data[i].nO, params, energy, grad);
			cout << "Energy: " << energy << endl;
		}
	}
	cout << "Time: " << (getTime() - startT)/1000.0 << "ms"  << endl;
	return 0;
	
}

#include<iostream.h>
#include<stdlib.h>
#include<cstring>
#include<fstream.h>
#include<ga/ga.h>

int main() {
	float mach, aoa;
	int selection = 1;
	int i,j;
	float uncert_min, uncert_max, float resultArray[1000], worst = 0, best = 10;

	if(argc>1) selection=atoi(argv[1]);

	char argu[40], argu2[20];
	switch(selection) {
		//only varying mach
		case 1:
			ofstream varyMach("varymach.dat", ios::trunc);
			varyMach.close();
		//varying both mach and aoa

			aoa = 2.0;

			ofstream varyMach("varymach.dat", ios::app);
			for(i=1; i<= 10;i++) {
				uncert_min = 0.5-0.5*(i*0.1);
				uncert_max = 0.5+0.5*(i*0.1);


				varyMach<<"Varying mach from "<<uncert_min<<" to "<<uncert_max<<" :"<<endl;
				for(j=0; j<1000; j++) {
					argu = "./lana ";
					mach = GARandomFloat(uncert_min, uncert_max);
					sprintf(argu2,"%f %f", mach, aoa);
					strcat(argu, argu2);
					system(argu);
					ifstream resultFile("result.dat");
					resultFile>>resultArray[j];
					if(resultArray[j]<best) best = resultArray[j];
					if(resultArray[j]>worst) worst = resultArray[j];
					varyMach<<j<<". "<<resultArray[j]<<endl;
				}
				varyMach<<"Best : "<<best<<endl;
				varyMach<<"Worst: "<<worst<<endl;
			}

			varyMach.close();


			break;

		//only varying aoa
		case 2: break;

		//varying both mach and aoa
		case 3: break;

		default: break;


	}
	return 0;
}
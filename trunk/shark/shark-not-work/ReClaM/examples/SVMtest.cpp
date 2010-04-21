
#include <Rng/GlobalRng.h>
#include <ReClaM/svm.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <iostream>


using namespace std;


int main()
{
	unsigned int e;
	double a, v, rad;
	unsigned int s;
	Rng::seed(10);

	double C = 5.0;
	double sigma = 1.0;
	double gamma = 0.5 / (sigma * sigma);
	RBFKernel k(gamma);

	cout << endl;
	cout << "We use a radial bases function kernel with a width of sigma=" << sigma << endl;
	cout << "and a complexity penalty of C=" << C << " for the support vector machine." << endl;
	cout << endl;

	// create the xor problem with uniformly distributed examples
	unsigned int examples = 100;
	cout << "Please enter a number of examples (20 to 200): " << flush;
	scanf("%d", &examples);
	getchar();		// read the '\n' character
	if (examples < 20) examples = 20;
	if (examples > 200) examples = 200;

	Array<double> x(examples, 2);
	Array<double> y(examples);
	cout << "Generating " << examples << " examples from the XOR distribution ..." << flush;
	for (e=0; e<examples; e++)
	{
		s = 0;
		a = Rng::uni(-2.5, 2.5);
		x(e, 0) = a;
		if (a >= 0.0) s++;
		a = Rng::uni(-1.75, 1.75);
		x(e, 1) = a;
		if (a >= 0.0) s++;
		if (s & 1) y(e) = -1.0; else y(e) = +1.0;
	}
	cout << " done." << endl;

	// create the SVM for prediction
	SVM svm(&k, x, y, false);

	// create a training scheme and an optimizer for learning
	C_SVM Csvm(&svm, C, C, true);
	SVM_Optimizer SVMopt;
	SVMopt.init(Csvm);

	// train the SVM
	cout << "Support Vector Machine training ..." << flush;
	SVMopt.optimize(svm, SVM_Optimizer::dummyError, x, y);
	cout << " done." << endl;

	// output the non-trivial components of alpha
	cout << endl << "The support vector indices are given in square brackets, each" << endl;
	cout << "followed by the value of the corresponding alpha coefficient:" << endl;
	unsigned int u, uc = x.dim(0);
	unsigned int SV = 0;
	unsigned int BSV = 0;
	for (u=0; u<uc; u++)
	{
		v = svm.get_alpha(u);
		if (v > 0.0)
		{
			SV++;
			printf("[%3.3d] ->%7.5g ", u, v);
			if (v >= C) BSV++;
		}
	}
	printf("\nSolution offset: b = %8.6g\n", svm.get_b());
	cout << "Number of support vectors: " << SV << endl;
	cout << "Number of bounded support vectors: " << BSV << endl;
	cout << endl;

	// Let the SVM predict a grid of test data to
	// estimate its accuracy on the XOR distribution
	int xx, yy;
	int acc = 0;
	Array<double> test(2);
	Array<double> output;
	for (yy=0; yy<100; yy++)
	{
		// numbers from [-1.75, 1.75]
		test(1) = (yy - 49.5) * 0.035;
		for (xx=0; xx<100; xx++)
		{
			// numbers from [-2.5, 2.5]
			test(0) = (xx - 49.5) * 0.05;
			svm.model(test, output);
			if (output(0) * test(0) * test(1) > 0.0) acc++;
		}
	}
	cout << "Estimated accuracy: " << 0.01 * acc << "%" << endl << endl;

	// output the solution to gnuplot
	char gp[65536]; gp[0] = 0;
	strcat(gp, "set title \"SVM example\"\n");
	strcat(gp, "unset key\n");
	strcat(gp, "unset label\n");
	strcat(gp, "unset clabel\n");
	strcat(gp, "set rmargin 0\n");
	strcat(gp, "set lmargin 0\n");
	strcat(gp, "set tmargin 0\n");
	strcat(gp, "set bmargin 0\n");
	strcat(gp, "set multiplot\n");
	strcat(gp, "set origin 0.0, 0.0\n");
	strcat(gp, "set size 1.0, 1.0\n");
	strcat(gp, "set palette gray\n");
	strcat(gp, "set pm3d map\n");
	strcat(gp, "set xrange [-2.5:2.5]\n");
	strcat(gp, "set yrange [-1.75:1.75]\n");
	strcat(gp, "set zrange [-10:10]\n");
	strcat(gp, "set samples 100\n");
	strcat(gp, "set isosamples 100\n");
	strcat(gp, "f(x,y)=");
	char part[256];
	bool first = true;
	for (u=0; u<uc; u++)
	{
		v = svm.get_alpha(u);
		if (v > 0.0)
		{
			if (! first && y(u)*v > 0.0) strcat(gp, "+");
			first = false;
			sprintf(part, "%g*exp(-%g*((x-%g)**2+(y-%g)**2))", y(u)*v, gamma, x(u, 0), x(u, 1));
			strcat(gp, part);
		}
	}
	if (svm.get_b() >= 0.0) strcat(gp, "+");
	sprintf(part, "%g\n", svm.get_b());
	strcat(gp, part);
	strcat(gp, "splot f(x,y)\n");
	strcat(gp, "unset pm3d\n");
	strcat(gp, "set contour\n");
	strcat(gp, "unset surface\n");
	strcat(gp, "set cntrparam levels incr 0,1,0\n");
	strcat(gp, "splot f(x,y)\n");
	strcat(gp, "set cntrparam levels incr -1,1,-1\n");
	strcat(gp, "splot f(x,y) with dots\n");
	strcat(gp, "set cntrparam levels incr 1,1,1\n");
	strcat(gp, "splot f(x,y) with dots\n");
	strcat(gp, "unset contour\n");
	strcat(gp, "set surface\n");
	strcat(gp, "set parametric\n");
	strcat(gp, "set isosamples 25\n");
	for (u=0; u<uc; u++)
	{
		v = svm.get_alpha(u);
		if (v > 0.0)
		{
			rad = 0.05 * sqrt(v / C);
			sprintf(part, "x(u,v)=%g+%g*cos(u)*cos(v)\ny(u,v)=%g+%g*sin(u)*cos(v)\nz(u,v)=%g*sin(v)\nsplot x(u,v),y(u,v),z(u,v)\n", x(u, 0), rad, x(u, 1), rad, rad);
			strcat(gp, part);
		}
	}
	strcat(gp, "splot u,0,v\n");
	strcat(gp, "splot 0,u,v\n");
	strcat(gp, "unset parametric\n");

	FILE* pipe = popen("gnuplot", "w");
	if (pipe != NULL)
	{
		fwrite(gp, 1, strlen(gp), pipe);
		fflush(pipe);

		cout << "The gnuplot window shows the SVM solution as a gray-coded" << endl;
		cout << "function on the input space. This function is a non-linear" << endl;
		cout << "pull-back of the linear SVM solution in a reproducing" << endl;
		cout << "kernel Hilbert space implicitly defined by the kernel." << endl;
		cout << "The cross shows the XOR problem class boundaries, while" << endl;
		cout << "the solid curved line shows the SVM decision boundary. The" << endl;
		cout << "dotted curves are the +1 and the -1 niveaus on which all" << endl;
		cout << "unbounded support vectors are located." << endl;
		cout << "The dots represent the support vectors. Their surface" << endl;
		cout << "areas represent the sizes of the corresponding solution" << endl;
		cout << "coefficients." << endl;
		cout << endl;
		cout << "*** press enter to quit ***" << flush;
		getchar();

		wait4(-2, NULL, 0, NULL);
		pclose(pipe);
	}
	else
	{
		cout << "*** unable to call gnuplot ***" << endl;
	}

	return 0;
}

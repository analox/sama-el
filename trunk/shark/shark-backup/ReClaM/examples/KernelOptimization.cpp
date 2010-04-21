
#include <Rng/GlobalRng.h>
#include <ReClaM/svm.h>
#include <ReClaM/KTA.h>
#include <ReClaM/RadiusMargin.h>
#include <ReClaM/optRprop.h>


int main()
{
	try{

	// create the chess board example
	unsigned int examples = 1000;
	unsigned int dim = 2;
	double C = 1000.0;
	double sigma = 0.5;

	Array<double> x(examples, dim);
	Array<double> y(examples);

	unsigned int e, d;
	double a;
	unsigned int s;
	Rng::seed(10);
	for (e=0; e<examples; e++)
	{
		s = 0;
		for (d=0; d<dim; d++)
		{
			a = Rng::uni(0.0, 4.0);
			x(e, d) = a;
			s += (int)a;
		}
		if (s & 1) y(e) = -1.0; else y(e) = +1.0;
	}

	// create two svm models based on kernel models
	RBFKernel* k[2];
	k[0] = new RBFKernel(0.5 / (sigma * sigma));
	k[1] = new RBFKernel(0.5 / (sigma * sigma));
	Model* model[2];
	SVM svm(k[0], x, y);
	model[0] = new C_SVM(&svm, C, C);
	model[1] = k[1];

	// create two error functions
	ErrorFunction* err[2];
	err[0] = new RadiusMargin();
	err[1] = new negativeKTA();

	// create an optimizer
	IRpropPlus* optimizer[2];
	optimizer[0] = new IRpropPlus();
	optimizer[1] = new IRpropPlus();

	// define names for the output
	char name[2][30] = {"Radius Margin Quotient", "Kernel Target Alignment"};

	// do two independent optimization runs
	int j, i;
	for (j=0; j<2; j++)
	{
		cout << endl << "Optimization run " << (j+1) << ":" << endl;
		cout << "Kernel optimization using the " << name[j] << "." << endl;
		cout << "initial kernel parameter: gamma = " << k[j]->getParameter(0) << endl;
		cout << "  Initializing the optimizer ..." << flush;
		optimizer[j]->init(*model[j]);
		cout << " done." << endl;

		for (i=0; i<30; i++)
		{
			cout << "  optimization iteration " << (i+1) << flush;
			double f = optimizer[j]->optimize(*model[j], *err[j], x, y);
			cout << "   f = " << f << endl;
		}
		cout << "final kernel parameter: gamma = " << k[j]->getParameter(0) << endl;
	}

	// clean up
	delete k[0];
	delete k[1];
	delete err[0];
	delete err[1];
	delete model[0];
	delete optimizer[0];
	delete optimizer[1];

	}
	catch (const char* e)
	{
		cout << "EXCEPTION: " << e << endl;
	}

	return 0;
}

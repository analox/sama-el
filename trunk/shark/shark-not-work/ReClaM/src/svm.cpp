
#include <math.h>
#include <ReClaM/svm.h>
#include <ReClaM/C_Solver.h>
#include <fstream>
#include <iostream>
// #include <iomanip>


using namespace std;


////////////////////////////////////////////////////////////////////////////////


SVM::SVM(KernelFunction* pKernel, bool bSignOutput)
{
	kernel = pKernel;
	signOutput = bSignOutput;

	x = NULL;
	y = NULL;
	bOwnMemory = false;

	examples = 0;
	dimension = 0;
}

SVM::SVM(KernelFunction* pKernel, const Array<double>& input, const Array<double>& target, bool bSignOutput)
{
	kernel = pKernel;
	signOutput = bSignOutput;

	setTrainingData(input, target);
}

SVM::~SVM()
{
	if (bOwnMemory)
	{
		delete x;
		delete y;
	}
}


void SVM::setTrainingData(const Array<double>& input, const Array<double>& target)
{
	x = &input;
	y = &target;
	bOwnMemory = false;

	examples = input.dim(0);
	dimension = input.dim(1);
	SIZE_CHECK(target.dim(0) == examples);

	parameter.resize(examples + 1);
	parameter = 0.0;
}

void SVM::model(const Array<double>& input, Array<double>& output)
{
	const double* xx;
	const double* ii;
	unsigned int i;
	double a;
	double alpha;
	if (input.ndim() == 1)
	{
		output.resize(1, false);
		a = parameter(examples);
		xx = x->elemvec();
		ii = input.elemvec();
		for (i=0; i<examples; i++)
		{
			alpha = parameter(i);
			if (alpha > 0.0)
			{
				a += alpha * (*y)(i) * kernel->eval(xx, ii, dimension);
			}
			xx += dimension;
		}
		output(0) = a;
	}
	else if (input.ndim() == 2)
	{
		int j, jc = input.dim(0);
		output.resize(jc, false);
		for (j=0; j<jc; j++)
		{
			a = parameter(examples);
			xx = x->elemvec();
			ii = input[j].elemvec();
			for (i=0; i<examples; i++)
			{
				alpha = parameter(i);
				if (alpha > 0.0)
				{
					a += alpha * (*y)(i) * kernel->eval(xx, ii, dimension);
				}
				xx += dimension;
			}
			output(j) = a;
		}
	}
	else throw "[SVM::model] invalid dimension";
}

void SVM::modelDerivative(const Array<double>& input, Array<double>& derivative)
{
	const double* xx;
	const double* ii;
	unsigned int i;
	double v;
	if (input.ndim() == 1)
	{
		derivative.resize(1, examples + 1, false);
		derivative = 0.0;
		xx = x->elemvec();
		ii = input.elemvec();
		for (i=0; i<examples; i++)
		{
			v = (*y)(i) * kernel->eval(xx, ii, dimension);
			xx += dimension;
			derivative(0, i) = v;
		}
		derivative(0, examples) = 1.0;
	}
	else if (input.ndim() == 2)
	{
		int j, jc = input.dim(0);
		derivative.resize(jc, examples + 1, false);
		derivative = 0.0;
		for (j=0; j<jc; j++)
		{
			xx = x->elemvec();
			ii = input[j].elemvec();
			for (i=0; i<examples; i++)
			{
				v = (*y)(i) * kernel->eval(xx, ii, dimension);
				xx += dimension;
				derivative(j, i) = v;
			}
			derivative(j, examples) = 1.0;
		}
	}
	else throw "[SVM::modelDerivative] invalid dimension";
}

void SVM::modelDerivative(const Array<double>& input, Array<double>& output, Array<double>& derivative)
{
	const double* xx;
	const double* ii;
	unsigned int i;
	double a, v;
	if (input.ndim() == 1)
	{
		derivative.resize(1, examples + 1, false);
		derivative = 0.0;
		output.resize(1, false);
		a = parameter(examples);
		xx = x->elemvec();
		ii = input.elemvec();
		for (i=0; i<examples; i++)
		{
			v = (*y)(i) * kernel->eval(xx, ii, dimension);
			xx += dimension;
			a += parameter(i) * v;
			derivative(0, i) = v;
		}
		output(0) = a;
		derivative(0, examples) = 1.0;
	}
	else if (input.ndim() == 2)
	{
		int j, jc = input.dim(0);
		output.resize(jc, false);
		derivative.resize(jc, examples + 1, false);
		derivative = 0.0;
		for (j=0; j<jc; j++)
		{
			a = parameter(examples);
			xx = x->elemvec();
			ii = input[j].elemvec();
			for (i=0; i<examples; i++)
			{
				v = (*y)(i) * kernel->eval(xx, ii, dimension);
				xx += dimension;
				a += parameter(i) * v;
				derivative(j, i) = v;
			}
			output(j) = a;
			derivative(j, examples) = 1.0;
		}
	}
	else throw "[SVM::modelDerivative] invalid dimension";
}

bool SVM::LoadSVMModel(std::istream& is)
{
	char buffer[50];
	unsigned int t, d;

	// read the header line
	is.read(buffer, 21);

	buffer[21] = 0;
	if (strcmp(buffer, "Shark SVM model\r\nSV: ") != 0) return false;
	if (! is.good()) return false;

	// read the number of support vectors
	is >> examples;

	is.read(buffer, 7); buffer[7] = 0;
	if (strcmp(buffer, "\r\nDIM: ") != 0) return false;
	if (! is.good()) return false;

	// read the input space dimension
	is >> dimension;

	is.read(buffer, 10); buffer[10] = 0;
	if (strcmp(buffer, "\r\nkernel: ") != 0) return false;
	if (! is.good()) return false;

	// read the kernel parameters
	is >> (*kernel);

	is.read(buffer, 14); buffer[14] = 0;
	if (strcmp(buffer, "coefficients: ") != 0) return false;
	if (! is.good()) return false;

	// read alpha and b
	parameter.resize(examples + 1);
	for (t=0; t<examples; t++)
	{
		is >> parameter(t);
		is.read(buffer, 1);
		if (buffer[0] != ' ') return false;
	}
	is >> parameter(examples);
	is.read(buffer, 2); buffer[2] = 0;
	if (strcmp(buffer, "\r\n") != 0) return false;
	if (! is.good()) return false;

	// read the support vectors with labels
	Array<double>* sv = new Array<double>(examples, dimension);
	Array<double>* l = new Array<double>(examples);
	for (t=0; t<examples; t++)
	{
		for (d=0; d<dimension; d++)
		{
			is >> (*sv)(t, d);
			is.read(buffer, 1);
			if (buffer[0] != ' ') return false;
		}
		is >> (*l)(t);
		is.read(buffer, 2); buffer[2] = 0;
		if (strcmp(buffer, "\r\n") != 0) return false;
		if (! is.good()) return false;
	}

	bOwnMemory = true;
	x = sv;
	y = l;

	return true;
}

bool SVM::SaveSVMModel(std::ostream& os)
{
	unsigned int t, d;
	unsigned int T = 0;
	for (t=0; t<examples; t++) if (parameter(t) > 0.0) T++;

	// write the header line
	os.write("Shark SVM model\r\nSV: ", 21);

	// write the number of support vectors
	os << T;

	os.write("\r\nDIM: ", 7);
	if (! os.good()) return false;

	// write the input space dimension
	os << dimension;

	os.write("\r\nkernel: ", 10);
	if (! os.good()) return false;

	// write the kernel parameters
	os << (*kernel);

	os.write("coefficients: ", 14);
	if (! os.good()) return false;

	// write alpha and b
	for (t=0; t<examples; t++) if (parameter(t) > 0.0) os << parameter(t) << " ";
	os << parameter(examples) << "\r\n";
	if (! os.good()) return false;

	// write the support vectors with labels
	for (t=0; t<examples; t++)
	{
		if (parameter(t) > 0.0)
		{
			for (d=0; d<dimension; d++)
			{
				os << (*x)(t, d) << " ";
			}
			os << (*y)(t);
			os.write("\r\n", 2);
			if (! os.good()) return false;
		}
	}

	return true;
}


////////////////////////////////////////////////////////////////////////////////


C_SVM::C_SVM(SVM* pSVM, double Cplus, double Cminus, bool norm2, bool unconst)
{
	svm = pSVM;
	C_plus = Cplus;
	C_minus = Cminus;
	C_ratio = Cminus / Cplus;
	norm2penalty = norm2;
	exponential = unconst;

	KernelFunction* kernel = svm->getKernel();
	unsigned int k, kc = kernel->getParameterDimension();
	parameter.resize(kc + 1, false);
	parameter(0) = (exponential) ? log(Cplus) : Cplus;
	for (k=0; k<kc; k++) parameter(k + 1) = kernel->getParameter(k);
}

C_SVM::~C_SVM()
{
}

void C_SVM::model(const Array<double>& input, Array<double> &output)
{
	svm->model(input, output);
}

void C_SVM::setParameter(unsigned int index, double value)
{
	Model::setParameter(index, value);
	if (index == 0)
	{
		// set C_+ and C_-
		if (exponential)
		{
			value = exp(value);
		}
		C_plus = value;
		C_minus = C_ratio * value;
	}
	else
	{
		// set a kernel parameter
		svm->getKernel()->setParameter(index - 1, value);
	}
}

double C_SVM::get_Cplus()
{
	return C_plus;
}

double C_SVM::get_Cminus()
{
	return C_minus;
}


////////////////////////////////////////////////////////////////////////////////


SVM_Optimizer::SVM_Optimizer()
{
	solver = NULL;
	printInfo = false;
	cacheMB = 100;
}

SVM_Optimizer::~SVM_Optimizer()
{
	if (solver != NULL) delete solver;
}


void SVM_Optimizer::init(Model& model)
{
	C_SVM* pSVM = dynamic_cast<C_SVM*>(&model);
	if (pSVM == NULL)
		throw "[SVM_Optimizer::init] The model is not a valid C_SVM.";

	Cplus = pSVM->get_Cplus();
	Cminus = pSVM->get_Cminus();
	norm2penalty = pSVM->is2norm();
}

void SVM_Optimizer::init(C_SVM& model, bool verbose, unsigned int mb)
{
	init(model);
	printInfo = verbose;
	cacheMB = mb;
}

double SVM_Optimizer::optimize(Model& model, ErrorFunction& error, const Array<double>& input, const Array<double>& target)
{
	SVM* pSVM = dynamic_cast<SVM*>(&model);
	if (pSVM == NULL)
	{
		C_SVM* pC_SVM = dynamic_cast<C_SVM*>(&model);
		if (pC_SVM != NULL) pSVM = pC_SVM->getSVM();
	}

	if (pSVM == NULL)
		throw "[SVM_Optimizer::optimize] The model is not a valid SVM or C_SVM.";

	optimize(*pSVM, input, target);

	return 0.0;
}

void SVM_Optimizer::optimize(SVM& model, const Array<double>& input, const Array<double>& target)
{
	if (solver != NULL) delete solver;

	model.setTrainingData(input, target);

	KernelFunction* kernel = model.getKernel();
	unsigned int i, examples = model.get_examples();

	Array<double> ones(examples);
	Array<double> C(examples);
	ones = 1.0;
	for (unsigned int e=0; e<examples; e++) C(e) = (target(e) > 0.0) ? Cplus : Cminus;

	Array<double> alpha(examples);
	for (i=0; i<examples; i++) alpha(i) = model.getParameter(i);
	double b = model.getParameter(examples);

	solver = new C_Solver(kernel, input, target, ones, target, C, alpha, b, printInfo, cacheMB, 0.001, "", true, norm2penalty);

	for (i=0; i<examples; i++) model.setParameter(i, alpha(i));
	model.setParameter(examples, b);
}


// static dummy member
ErrorFunction& SVM_Optimizer::dummyError = *((ErrorFunction*)NULL);

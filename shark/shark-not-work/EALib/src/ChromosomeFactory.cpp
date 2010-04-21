
#include <EALib/ChromosomeFactory.h>
#include <EALib/ChromosomeT.h>


Chromosome* CreateChromosome(const char* type)
{
	if (strcmp(type, typeid(double).name()) == 0)
	{
		return new ChromosomeT<double>();
	}
	else if (strcmp(type, typeid(int).name()) == 0)
	{
		return new ChromosomeT<int>();
	}
	else if (strcmp(type, typeid(bool).name()) == 0)
	{
		return new ChromosomeT<bool>();
	}
	else if (strcmp(type, typeid(char).name()) == 0)
	{
		return new ChromosomeT<char>();
	}

	printf("\n\n[CreateChromosome] unknown type: %s\n", type);
	exit(1);
	return NULL;
}

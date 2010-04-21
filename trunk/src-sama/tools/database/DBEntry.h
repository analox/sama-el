
#ifndef DBEntry_H
#define DBEntry_H
#include <EALib/Individual.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

class DBEntry {
public:
	int t;
	Individual startPoint, endPoint;
	double cost;
	double fi;
	DBEntry(int t, Individual s, Individual e, double c, double f);
	DBEntry(const DBEntry & entry);
	DBEntry();
	int getGenNum() {return t;}
	double getCost() {return cost;}
	Individual & getStartPoint() {return startPoint;}
	Individual & getEndPoint() {return endPoint;}
	string toString();
};
#endif

#include "DBEntry.h"

	DBEntry::DBEntry(int t, Individual s, Individual e, double c, double f) {
		this->t = t;
		this->startPoint = s;
		this->endPoint = e;
		this->cost = c;
		this->fi= f;
	}

	DBEntry::DBEntry() {}

	DBEntry::DBEntry(const DBEntry & entry)
	{
		this->t = entry.t;
		this->startPoint = entry.startPoint;
		this->endPoint = entry.endPoint;
		this->cost = entry.cost;
		this->fi = entry.fi;
	}

	string DBEntry::toString() {
		ostringstream res;
		res << "=========== ";
		res << t << endl;
		res << startPoint[0] << endl;
		res << endPoint[0] << endl;
		res << "FI: " << fi << ", Cost: " << cost << endl;
		return res.str();
	}

	// UNIT TEST
	int _DBEntry_main()
	{
		cout << "Testing DBEntry" << endl;
		double _a[] = {1, 2, 3};
		vector<double> a(_a, _a + 3);
		//for (int i = 0; i < a.size(); i++) 
		//	cout << a[i] << (i == (a.size()-1 ) ? " " : ", ");
		//cout << endl;

		ChromosomeT<double> chrom1(a);
		Individual startPoint(chrom1);
		//cout << startPoint[0] << endl;

		double _b[] = {4, 6, 3};
		vector<double> b(_b, _b + 3);
		ChromosomeT<double> chrom2(b);
		Individual endPoint(chrom2);
		//cout << endPoint[0] << endl;

		DBEntry entry1(0, startPoint, endPoint, 10, 0.123);
		cout << entry1.toString();
		DBEntry entry2(entry1);
		entry2.t++;
		cout << entry2.toString();
		return 0;

	}


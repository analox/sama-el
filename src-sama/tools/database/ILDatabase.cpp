#include "ILDatabase.h"

/**
 * Class to manage list of (starting, end point) for different 
 * individual learning
 * @author Le Minh Nghia, NTU-Singapore
 *
 */	
	ILDatabase::ILDatabase() {}
	ILDatabase::ILDatabase(int n)
	{
		if (n > 0) 
		{
			nMethods = n;
			DB.resize(nMethods);;
			listArray.resize(nMethods);
			for (int i = 0; i < nMethods; i++) 
			{
				list<DBEntry> l;
				DB.push_back(l);
			}
		}
		else
			cout << "Invalid input..." << endl;
	}
	/**
	 * Get list of entries for methodID
	 * @param methodID
	 * @return
	*/ 
	list<DBEntry> & ILDatabase::getListForMethod(int methodID)
	{
		return DB[methodID];
	}
	
	/**
	 * Add a record of some methodID to list
	 * @param entri
	 * @param methodID
	 */
	void ILDatabase::add(DBEntry entri, int methodID)
	{
		list<DBEntry> & aList = this->getListForMethod(methodID);
		if (aList.size() >= MAX_ENTRIES_PER_LIST)
			aList.pop_front();
		aList.push_back(entri);
	}
	/*  */
	void ILDatabase::prepareListArray()
	{
		for (int i = 0; i <  nMethods; i++) {
			list<DBEntry> & listDB = this->getListForMethod(i);
			// check if list empty
			if (!listDB.empty()) 
			{
				// Get chromosome dimension
				list<DBEntry>::iterator lIter = listDB.begin();	
				int nDim = ((*lIter).getStartPoint())[0].size();

				_Matrix<double> mat(listDB.size(), nDim);

				// get array of start points
				int indexDB = 0;
				for (lIter=listDB.begin(); lIter != listDB.end(); ++lIter)
				{
					for (int dim = 0; dim < nDim; dim++)
						mat.at(indexDB, dim) = dynamic_cast<ChromosomeT<double>&>(((*lIter).getStartPoint())[0]).at(dim);
					indexDB++;
				}
				// add to list		
				listArray.at(i) = mat;	
			}
		}
	}

	/**
	 * Get list of start points for methodID
	 * @param methodID
	 * @return
	*/ 

	_Matrix<double> & ILDatabase::getListArrayForMethod(int methodID)
	{
		return listArray[methodID];
	}
	
	/**
	 * Reset the database
	*/ 
	void ILDatabase::resetDB()	{
		for (int i = 0; i < nMethods; i++) 
		{
			list<DBEntry> & l = DB[i];
			l.clear();
			listArray.clear();
		}
	}
	
	/* */
	string ILDatabase::toString()
	{
		ostringstream res;
		for (int i = 0; i < nMethods; i++)
		{
			res << "\t List ID = " << i << ", " << DB[i].size() << " entries \n";
			  
			//list<DBEntry>::iterator lIter;

                 	//for(lIter=DB[i].begin(); lIter != DB[i].end(); ++lIter) 
			//	res << (*lIter).toString() << endl;	
		}
		return res.str();
	}
	
	// UNIT TEST
	int _ILDatabase_main() 
	{
		cout << "ILDatabase test" << endl;
		ILDatabase db(2);

		double _a[] = {1, 2, 3};
		vector<double> a(_a, _a + 3);
		ChromosomeT<double> chrom1(a);
		Individual startPoint(chrom1);

		double _b[] = {4, 6, 3};
		vector<double> b(_b, _b + 3);
		ChromosomeT<double> chrom2(b);
		Individual endPoint(chrom2);

		DBEntry entry1(0, startPoint, endPoint, 10, 0.123);
		DBEntry entry2(entry1); entry2.t++;

		db.add(entry1, 0); db.add(entry2,1);
		cout << db.toString();

		db.prepareListArray();
		cout << "Array type:" << endl;
		cout << db.getListArrayForMethod(0);
		cout << "Array type:" << endl;
		cout << db.getListArrayForMethod(1);
		return 0;
	}


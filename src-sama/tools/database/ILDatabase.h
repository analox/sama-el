#ifndef ILDatabase_H
#define ILDatabase_H
/**
 * Class to manage list of (starting, end point) for different 
 * individual learning
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
#include "DBEntry.h"
#include <list>
#include <Array/ArrayTable.h>
#include "../Matrix.cpp"

using namespace std;

class ILDatabase {
public:
	static const int MAX_ENTRIES_PER_LIST = 200;
	int nMethods;
	vector<list<DBEntry> > DB;
	vector<_Matrix<double> > listArray;

	ILDatabase();	
	ILDatabase(int n);
	/**
	 * Get list of entries for methodID
	 * @param methodID
	 * @return
	 */
	list<DBEntry> & getListForMethod(int methodID);
	/**
	 * Add a record of some methodID to list
	 * @param entri
	 * @param methodID
	 */
	void add(DBEntry entri, int methodID);
	void prepareListArray();
	/**
	 * Get list of start points for methodID
	 * @param methodID
	 * @return
	 */
	_Matrix<double> & getListArrayForMethod(int methodID);
	/**
	 * Reset the database
	 */
	void resetDB();
	string toString();
};
#endif

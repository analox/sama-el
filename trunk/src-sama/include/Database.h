#ifndef __DATABASE__
#define __DATABASE__

#include "Array/Array.h"
#include <vector>
#include <string>

/*! \file */ 
//! This class provides the general functionalities for a simple database. 

/*!The database is based on two Arrays containing the Input and Target data sets. It stores data sets in a chronological order. The maximum number of data sets is restricted to a number of MaxTrainingData data sets, where each data set consists of a vector of input data and a vector of target data which can be used for several approximation models (e.g. neural network, regression models, kriging model). Since it is beneficial to use only a limited number of the newest training data (especially in case of high dimensional input data where a local approximation of the fitness function near the actual position of the population is needed) the number of data sets which are used for the training is limited to NumTrainingData. For the training only NumTrainingData data sets are used. The restriction to a limited number is mainly introduced due to stability reasons and the value of MaxTrainingData can usually set to sufficiently large value. NewestEntry gives the position of the last added data set. */

class Database {

 public:

 private:

  //! actual number of available training data 
  unsigned NumTrainingData;
  //! maximum number of training data
  unsigned MaxTrainingData;
  //! number of used training data
  unsigned UsedTrainingData;
  //! maximum length of training data archive
  int MaxArchiveLength;
  //! position of the newest dataset in archive
  int NewestEntry;

  //! array containing input data of given data sets
  Array < double > Input;
  //! array containing target data of given data sets
  Array < double > Target;
  //! array conotaining informations on the generation
  Array < int > Information;

  // ======================================================================
  
 public:
  
  /*! \brief constructor */
  Database();

  /*! \brief This function adds a new set of training data to the database.
      \param NewInput Array containing new Input data
      \param NewTarget Array containing new Target data
      \retval 0 ok
      \retval 1 if dimensions of the new data set do not fit to an existing database
  */
  int AddTrainingData( Array < double > NewInput, 
		       Array < double > NewTarget );

  /*! \brief This function adds a new set of training data to the database with respect to the number of the generation.
      \param NewInput Array containing new Input data
      \param NewTarget Array containing new Target data
      \param numberOfGeneration number of the generation the data comes from
      \retval 0 ok
      \retval 1 if dimensions of the new data set do not fit to an existing database
  */
  int AddTrainingData ( Array < double > NewInput,
			Array < double > NewTarget,
			int numberOfGeneration );

  /*! \brief This function adds a new set of training data to the database with respect to the numbers of the generations.
      \param NewInput Array containing new Input data
      \param NewTarget Array containing new Target data
      \param Generation Array containing the number of the generation each data comes from
      \retval 0 ok
      \retval 1 if dimensions of the new data set do not fit to an existing database
  */
    int AddTrainingData ( Array < double > NewInput,
			  Array < double > NewTarget,
			  Array < int > Generation );
 
  /*! \brief This function adds one training data to the database.
      \param NewInput Vector containing new Input data
      \param NewTarget Vector containing new Target data
      \retval 0 ok
      \retval 1 if dimensions of the new data set do not fit to an existing database
  */
  int AddTrainingData ( std::vector < double > NewInput, 
			std::vector < double > NewTarget );

 
  /*! \brief This function adds one training data to the database with respect to the number of generation the data comes from.
      \param NewInput Vector containing new Input data
      \param NewTarget Vector containing new Target data
      \param numberOfGeneration number of the generation the data comes from
      \retval 0 ok
      \retval 1 if dimensions of the new data set do not fit to an existing database
  */
  int AddTrainingData ( std::vector < double > NewInput, 
			std::vector < double > NewTarget,
			int numberOfGeneration );

  /*! \brief This function deletes the first training data from the database.
      \retval 0 ok
      \retval 1 if database is empty
  */
  int DeleteFirstData();

  /*! \brief This function deletes the last training data from the database.
      \retval 0 ok
      \retval 1 if database is empty
  */
  int DeleteLastData();

  /*! \brief This function deletes training data at a given Index from the database.
      \param Index position of the data in the database which shall be deleted
      \retval 0 ok
      \retval 1 if database is empty
  */
  int DeleteDataAtPosition( int Index );

  /*! \brief This function deletes training data from a given Index1 to Index2 from the database.
      \param Index1 position of the data in the database from which shall be deleted
      \param Index2 position of the data in the database up to which shall be deleted
      \retval 0 ok
      \retval 1 if database is empty
  */
  int DeleteDataFromIndex1ToIndex2( int Index1,
				    int Index2 );

  /*! \brief This function deletes training data of generation <numberOfGeneration> from the database.
      \param numberOfGeneration number of the generation to be deleted
      \retval 0 ok
  */
  int DeleteGeneration( int numberOfGeneration );

  /*! \brief This function gets the last <NumTrainingData> data sets and stores them in the Arrays InputData and TargetData.
    
      \param InputData Array which contains the last <NumTrainingData> input data
      \param TargetData Array which contains the last <NumTrainingData> target data
      \retval 0 ok
      \retval 1 if database is empty
  */
  int GetLastData ( Array < double > &InputData, 
		    Array < double > &TargetData );

  /*! \brief This function gets the last <NumTrainingData> data sets and stores them in the Arrays InputData, TargetData and numberOfGeneration.
    
      \param InputData Array which contains the last <NumTrainingData> input data
      \param TargetData Array which contains the last <NumTrainingData> target data
      \param numberOfGeneration Array which contains the last <NumTrainingData> generation numbers of each data
      \retval 0 ok
      \retval 1 if database is empty
  */
  int GetLastData ( Array < double > &InputData, 
		    Array < double > &TargetData,
		    Array < int > &numberOfGeneration );

  /*! \brief This function gets the data at Index from the database and stores them in the Arrays InputData and TargetData.
    
      \param InputData Array which contains the input data at position Index
      \param TargetData Array which contains the target data at position Index
      \param Index position of data
      \retval 0 ok
      \retval 1 if Index is out of database
  */
  int GetDataAtPosition( Array < double > &InputData, 
			 Array < double > &TargetData,
			 int Index );

  /*! \brief This function gets the data at Index from the database and stores them in the Arrays InputData, TargetData and numberOfGeneration.
    
      \param InputData Array which contains the input data at position Index
      \param TargetData Array which contains the target data at position Index
      \param numberOfGeneration number of the generation the data comes from
      \param Index position of data
      \retval 0 ok
      \retval 1 if Index is out of database
  */
  int GetDataAtPosition( Array < double > &InputData, 
			 Array < double > &TargetData,
			 int &numberOfGeneration, 
			 int Index );

  /*! \brief This function gets the data for the generation <numberOfGeneration> from the database and stores them in the Arrays InputData and TargetData.
      \param numberOfGeneration number of the generation the data comes from
      \param InputData Array which contains the input data of generation <numberOfGeneration>
      \param TargetData Array which contains the target data of generation <numberOfGeneration>
      \param Index position of data
      \retval 0 ok
      \retval 1 if no data is found
  */
  int GetDataOfGeneration ( int numberOfGeneration,
			    Array < double >& InputData,
			    Array < double >& TargetData );

  /*! \brief This function loads training data from a file named Filename. If add is chosen to 0 consisting data in the database will be overwritten. If add equals 1 the loaded data will be appended to the database.
    
      \param Filename name of a file containing data sets
      \param add if add equals 0 the database is overwritten, if add equals 1 the data will be appended.
      \retval 0 ok
      \retval 1 if file cannot be opened
      \retval 2 if dimension of the loaded data does not fit to the existing database
  */
  int LoadTrainingData ( std::string Filename, 
			 int add);

  /*! \brief This function saves training data to a file named Filename.
      \param Filename name of a file in which the data sets are saved
      \retval 0 ok
      \retval 1 if file cannot be opened
      \retval 2 if database is empty
  */
  int SaveTrainingData ( std::string Filename );

  /*! \brief This function prints the data stored in the database on the screen.
      \retval 0 ok
  */
  int PrintTrainingData();

  /*! \brief This function returns the number of data in the database.
      \retval NumTrainingData number of data in the database
  */
  int GetNumTrainingData();

  /*! \brief This function returns the maximum number of data which can be used for training.
      \retval MaxTrainingData maximum number of data which may be used for training
  */
  int GetMaxTrainingData();

  /*! \brief This function returns the number of used training data.
      \retval UsedTrainingData number of data which is used for training
  */
  int GetUsedTrainingData();

  /*! \brief This function returns the maximum number of data which can be stored in the database.
      \retval MaxArchiveLength maximum number of data which can be stored in the database
  */
  int GetMaxArchiveLength();

  /*! \brief This function returns the position of the data which was recently added.
      \retval NewestEntry position of data which was added recently
  */
  int GetNewestEntry(); 

  /*! \brief This function sets the maximum number of data which can be stored in the database to NewLength.
      \param NewLength maximum number of data which can be stored
      \retval 0 ok
      \retval 1 if the new length is larger than the number of data in database. The oldest data will be removed until the length equals NewLength.
  */
  int SetMaxArchiveLength( int NewLength );
 
  /*! \brief This function sets the maximum number of training data to NewMax.
      \param NewMax maximum number of data which can be used for training
      \retval 0 ok
  */
  int SetMaxTrainingData( int NewMax );

  // ======================================================================

 private:

  int AddElement ( Array < double > NewInput, 
		   Array < double > NewTarget,
		   Array < int > Generation );

};

#endif

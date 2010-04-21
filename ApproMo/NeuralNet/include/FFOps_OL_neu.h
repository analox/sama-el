#ifndef FFOPS_OL_NEU_H
#define FFOPS_OL_NEU_H
#include "./FFOps_Basic_neu.h"
//#pragma interface
using namespace std;
//! This class is inherited from class \c FF_Basic_neu for mutations on a network
class FFOps_OL_neu : public FFOps_Basic_neu{
public:
  // modify connectivity
  //! Modify weights of a layered network
  /*!
    \param Ind representing the network to be modified
    \param prob1: probability for a connection between input and hidden nodes
    \param prob2: probability for a connection between hidden and output nodes
    \param prob3: probability for a connection between bias nodes and other nodes
*/
  void initConnections(Individual& Ind, double prob1=1.0, double prob2=1.0, double prob3=1.0);
  //! Add connections where possible and initialize the weight
  /*!
    \param Ind representing the network in which connection to be added
    \param low: lower bound of the interval for generating a random weight
    \param high: upper bound of the interval for generating a random weight
  */
  void addConnection(Individual& Ind, double low, double high);
  //! delete a connection 
  /*! \param Ind representing the network in which a connection is to be deleted 
      \param sigma ?
   */
  void deleteConnection(Individual& Ind, const double sigma=-1);
  //! delete a neuron
  /*! \param Ind representing the network in which a neuron is to be deleted 
   */
  void deleteNeuron(Individual& Ind);
  //! Add a neuron
  /*! Add connections between the newly generated neuron and input, output 
    neurons. Initialize the weights randomly with a given interval.
    \param Ind: individual representing the network to be processed
    \param low: lower bound of the interval for generating a random weight
    \param high: upper bound of the interval for generating a random weight
  */
  void addNeuron(Individual& Ind, double low, double high);
  //! Delete the output of a neuron
  void deleteNeuronGain(Individual& Ind, double sigma=-1);
  //! Delete the activation of a neuron
  void deleteNeuronActivation(Individual&, double=-1);
};
#endif

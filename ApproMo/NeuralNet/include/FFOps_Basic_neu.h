#ifndef FFOPS_BASIC_NEU_H
#define FFOPS_BASIC_NEU_H

#include "./Ops_Basic.h"
#include "./FFManip.h"

//#pragma interface
using namespace std;
//! Basic class for defining feed-forward neural networks
/*! This class defines feedforward neural networks
    based on the functions defined in \c Ops_basic.h and \c FFManip.h.
*/

class FFOps_Basic_neu : public Ops_Basic,
			protected FFManip{
public:
  /*__ Funktion kreiert ein neues Individuum nIn Eingangsneuronen
    (1. Parameter), nHid versteckten Neuronen (2. Parameter) und nOut
    Ausgangsneuronen (3. Parameter). Die Zahl der Zeitschichten wird
    somit auf 1 beschraenkt. */ 
  //! This function creats a generic neural network.
  /*! The neural network is generated given the number of input, 
    output and hidhen nodes.
      \param nIn: 
      number of input nodes
      \param nHidden: 
      number of hidden nodes
      \param nOut:
      number of output nodes.
  */
  Individual createEmpty(const unsigned nIn, const unsigned nHidden, 
			 const unsigned nOut);
  
  /*__ Funktion gibt die Verbindungsmatrix zurueck. Diese Matrix wird
    aber so organisiert, dass die Bias-Verbindungen sich in einer
    zusaetztlichen Spalte am Ende der Matrix befinden, also nicht mehr
    auf der Diagonalen. */
  //! This function returns the connection matrix of the neural network.
  /*! Notice that connections between the bias nodes
      are provided separately in column 
      at the end of the matrix.
      \param Ind: 
      Individual variable that defines the structure and
      parameters of the neural network.
  */  
  Array<int> getConnectionsReClaM(Individual& Ind);

  /*__ Funktion gibt die Verbindungsmatrix zurueck
    (vgl. getConnectionsReClaM), allerdings wird die zusaetzliche
    Spalte mit dem Bias nicht mit erzeugt, so dass sich eine
    quadratische Matrix ohne Bias-Information ergibt. */
  //! This function also returns connection matrix of the neural network.
  /*! However, the connection of the bias nodes are not included.
    \param Ind:
    Individual variable that defines the structure and
      parameters of the neural network.
  */
  Array<int> getConnectionsWithoutBias(Individual& Ind);

  /*__ Funktion gibt die Gewichtsmatrix zurueck. Diese Matrix wird
    aber so organisiert, dass die Gewichte der Bias-Verbindungen sich
    in einer zusaetztlichen Spalte am Ende der Matrix befinden, also
    nicht mehr auf der Diagonalen. */  
  //! This function returns the weigh value of the connections.
  /*! The weights connecting the bias nodes are listed
     in column at the end of the matrix.
   \param Ind: Individual variable that defines the structure and
      parameters of the neural network.
  */
  Array<double> getWeightsReClaM(Individual& Ind);

  /*__ Funktion setzt die Verbindungsmatrix des Individuums. Hierbei
    wird das ReClaMtypische Format verwendet, allerdings mit einer
    zusaetzlichen Bias-Spalte am Ende der Matrix. */ 
  //! This function defines the connections (structure) of a neural network.
  /*! The network structure is defined by the connection matrix.
    \param Ind: an individual representing a neural network
    \param cR: the connection matrix
  */

  void setConnectionsReClaM(Individual& Ind, Array<int> cR);

  /*__ */
  //! This function defines the weights of a neural network.
  /*! The network weights are given the connection matrix.
    \param Ind: an individual representing a neural network
    \param cW: the weight matrix
  */
  void setWeightsReClaM(Individual& Ind, Array<double> cW);
  //! Print out the connections as well as the weights of the network.
  /*! The data are printed out to a file whose name is supplied
    by \c filename.
    \param Ind: The network to be printed out, which is represented by
    an individual.
    \param filename: name of the file to store the connections and weights
     of the network.
  */ 
  void printNet(Individual& Ind, string filename);
  //! Print out the connections as well as the weights of the network.
  /*! The data are printed out to a file whose name is supplied
    by \c filename. However, the biases are not included.
    \param Ind: The network to be printed out, which is represented by
    an individual.
    \param filename: name of the file to store the connections and weights
     of the network.
  */
  void printNetWithoutBias(Individual& Ind, string filename);
  //! This function calculates the output of all neurons in the network.
  /*!
    \param Ind: The network to be calculated defined by an individual
  */
  Array<double> getGain(Individual& Ind);
  //! Set the activation value.
  /*!
    \param Ind: Individual representing the network
    \param activation: a matrix that contains the activation values
  */

  void setActivation(Individual& Ind, Array<double> activation);
  //! This function returns the activation values of a network
  /*! 
    \param Ind: Individual representing the network to be accessed
  */
  Array<double> getActivation(Individual& Ind);
};
#endif

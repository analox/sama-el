#ifndef OPS_BASIC_H
#define OPS_BASIC_H

#include <Array/Array.h>
#include <Array/ArrayIo.h>
#include <EALib/Individual.h>
#include <EALib/Chromosome.h>
#include <Rng/GlobalRng.h>
#include <Array/ArrayOp.h>
#include <stdio.h>
#include <string.h>
//#include <pvm3.h>
//#pragma interface
using namespace std;
//! This class supplies basic functions for getting access to neural netwoks stored in chromosome.
/*! Access operations include reading out the networks for mutation
    and again save the networks back to the chromosom.
*/
/*__ Diese Klasse stellt die wesentliche Funktionalitaet zur Verfuegung, 
     um Netze in die Chromosomen der Individuen der EALib zu speichern, aus
     diesen wieder auszulesen und Mutationen auszufuehren*/
class Ops_Basic{
 public:
  /*__ Kreiert ein neues Individuum. Die einzelnen Parameter geben die
       Zahl der Eingaenge nIn, Zahl der versteckten Neurone nHid, die
       Zahl der Ausgangsneurone nOut und die Zahl der
       unterschiedlichen Zeitschichten nDelay. Letztendlich werden
       zwei dreidimensionale Arrays mit den Dimensionen
       nDelay*(nIn*nHid*nOut)*(nIn*nHid*nOut) realisiert */
  //! Creat an new individual for representing a network.
  /*! Define a generic network by defining the number of
    input, hidden and output nodes. For generating recurrent networks,
    a parameter \n nDelay is also needed, whereas if feedforward networks
    are to be created, this should be fixed to 1.
    Finally, two three-dimensional matrices of the size \c nDelay*(nIn+nHid+nOut)*(nDelay+nIn+nHid+nOut) will be generated.
    
    \param nIn: number of input nodes
    \param nHid: number of hidden nodes
    \param nOut: number of output nodes
    \param nDelay: number of time delay
  */
  Individual createEmpty(const unsigned nIn, const unsigned nHid, const unsigned nOut, const unsigned nDelay);

  /*__ Funktion initialisiert alle Gewichte im Individuum. Die
    Gewichte werden gleichverteilt aus dem Intervall [low, high]
    gezogen, wobei low der zweite, high der dritte Parameter der
    Funktion ist */
  //! Initialize the weights randomly within a given interval.
  /*! 
    \param Ind: individual representing the network whose weights need
    to be initialized
    \param low: lower bound of the interval
    \param high: upper bound of the interval
  */

  void initWeights(Individual& Ind, double low, double high);

  /*__ Funktion gibt die Verbindungsmatrix in Form eines
    dreidimensionalen Arrays zurueck */
  //! Return the connection matrix of three dimensional. 
  /*
    \param Ind: individual representing the network
  */ 
  Array<int> getConnections(Individual& Ind); 
  
  /*__ Funktion gibt die Verbindungsmatrix nur einer Zeitschicht
    (zweiter Parameter) zurueck. Die Matrix ist quadratisch mit der
    Dimension (nIn*nHid*nOut). Die Zeitschichten werden von null
    beginned gezaehlt. */
  //! Returns the connection matrix of two dimensional of a given time
  /*!
    \param Ind: individual representing the network
    \param entry: time index for which connections should be returned.
  */

  Array<int> getConnections(Individual& Ind, unsigned entry);

  /*__ Funktion gibt die Gewichtsmatrix in Form eines
    dreidimensionalen Arrays zurueck */
  //! Return the weights of the network.
  /*!
     \param Ind: individual representing the network
  */
  Array<double> getWeights(Individual& Ind); 
  
  /*__ Funktion gibt die Gewichtsmatrix nur einer Zeitschicht (zweiter
    Parameter) zurueck. Die Matrix ist quadratisch mit der Dimension
    (nIn*nHid*nOut). Die Zeitschichten werden von null beginned
    gezaehlt. */
  //! Returns the weights of a given time
  /*!
    \param Ind: individual representing the network
    \param entry: time index for which weights should be returned.
  */
  Array<double> getWeights(Individual& Ind, unsigned entry);

  /*__ Funktion gibt die Zahl der Eingangsneurone des Individuums
    zurueck */
  //! Return the number of input neurons
  /*!
    \param Ind: individual representing the network to be accessed
  */
  unsigned getNIn(Individual& Ind);

  /*__ Funktion gibt die Zahl der versteckten Neurone des Individuums
    zurueck */
  //! Return the number of hidden neurons
  /*!
    \param Ind: individual representing the network to be accessed
  */
  unsigned getNHid(Individual& Ind);

  /*__ Funktion gibt die Zahl der Ausgangsneurone des Individuums
    zurueck */
  //! Return the number of output neurons
  /*!
    \param Ind: individual representing the network to be accessed
  */
  unsigned getNOut(Individual& Ind);

  /*__ Funktion gibt die Zahl der Zeitschichten des Individuums
     zurueck */
  //! Return the time index of the network
  /*!
    \param Ind: individual representing the network to be accessed
  */
  unsigned getNDelay(Individual& Ind);

  /*__ Funktion gibt die Zahl der freien Parameter, entsprechend der
     Anzahl der Eintraege in den Verbindungsmatrix, des Individuums
     zurueck */
  //! Return the number of all parameters contained in the network
  /*!
      \param Ind: individual representing the network to be accessed
  */
  unsigned getNParameter(Individual& Ind);

  /*__ Das Netz, welches im Individuum codiert ist, wir in den File
     gespeicher, welcher den Namen des zweiten Parameters traegt */
  //! Print the network structure and parameters to a file called \c filename
  /*
    \param Ind: individual representing the network to be printed
    \param filename: the file name in which the data will be printed
  */
  void printNet(Individual& Ind, string filename);

  /*__  Die Gewichte in der im letzten Parameter angegebenen
     Zeitschicht des Individuums wird durch die Werte der uebergebenen
     Matrix ersetzt. Die Zeitschichten werden von 0 beginnend
     gezaehlt. */
  //! Replace the weights of a given time index with those provided
  /*! 
    \param Ind: individual representing the network whose weights are to be
    replaced
    \param w: a weight matrix to replace the current weights
    \param entry: time index for the weights to be replaced
  */
  void setStructure(Individual& Ind, const Array<double>& w, unsigned entry);

  /*__  Die Gewichte des Individuums werden durch die Werte der
     uebergebenen Matrix ersetzt. */
  //! Replace all weights with those provided
  /*! 
    \param Ind: individual representing the network whose weights are to be
    replaced
    \param w: a weight matrix to replace the current weights
  */
  void setStructure(Individual& Ind, const Array<double>& w);

  /*__  Die Verbindungen in der im letzten Parameter angegebenen
     Zeitschicht des Individuums wird durch die Werte der uebergebenen
     Matrix ersetzt. Die Zeitschichten werden von 0 beginnend
     gezaehlt. */
  //! Replace the connections of a given time index with those provided
  /*! 
    \param Ind: individual representing the network whose connections 
    are to be replaced
    \param c: a connection matrix to replace the current connections
    \param entry: time index for the weights to be replaced
  */

  void setStructure(Individual& Ind, const Array<int>& c, unsigned entry);

  /*__  Die Verbindungen des Individuums werden durch die Werte der
     uebergebenen Matrix ersetzt. */
  //! Replace connections with those provided
  /*! 
    \param Ind: individual representing the network whose connections 
     are to be replaced
    \param c: a weight matrix to replace the current connections
  */
  void setStructure(Individual& Ind, const Array<int>& c);

  /*__ Die Gewichte werden normalverteilt mutiert. Der zweite
    Parameter gibt die Standardabweichung der Normalverteilung
    an. Dier 3. Parameter gibt an, welcher Anteil der Gewichte mutiert
    werden. */
  //! Mutate the weights
  /*! The weights are mutated using a normal distribution 
    with a given propabability.
    \param Ind: individual to be mutated
    \param stdv: standard diviation of the normal distribution
    \param prob:  mutation probability
  */
  void jogWeights(Individual& Ind, const double stdv, double prob= 1.0);

  /*__ Funktion transformiert alle Informationen des Individuums in
    eine PVM-taugliche Form. Dieser Befehl wirkt analog zu den
    sonstigen pvm_pk-Befehlen */
  //! Prepare a data structure ready for PVM
  /*!
     \param Ind: individual to be processed
  */
  // int pvm_pkIndividual(Individual& Ind);

  /*__ Funktion entpackt ein Individuum aus einem PVM-Stream */
  //! Receive data from a PVM stream 
  /*!
    \param Ind: individual to be received
  */
  //int pvm_upkIndividual(Individual& Ind);

protected:
  /*__ Funktion ermittelt die Gesamtzahl der Neurone eines
    Individuums */
  //! Return the sum of neuron number in the network 
    /*!
     \param Ind: individual to access
  */
  unsigned getNNeurons(Individual& Ind);
  
  /*__ Funktion veraendert die Zahl der versteckten Neurone auf den im
    zweiten Parameter angegebenen Wert. Es werden allerdings nur die
    Laengen der Chromosome veraendert, nicht aber die Gewichts- und
    Verbindungsmatrizen entsprechend umkopiert. */
  //! Change the number of hidden nodes according to a given number
  /*! Only the length of the chromosome is changes, the connection and weight
    matrices remain unchanged.
     \param Ind: individual to be changed
     \param nHid: intended number of hidden nodes
  */
  void setNHid(Individual& Ind, unsigned nHid);

  /*__ Ein verstecktes Neuron wird geloescht. Die versteckten Neurone
    werden, mit 0 beginnend durchgezaehlt, so dass der Wert 0 die
    Loeschung des 0. versteckten Neurons bewirkt. */
  //! Delete a hidden neuron according to a given number.
  /*! All hidden neurons are enumerated so that an arbitrary hidden neuron
    can be determined by an integer number.
    \param Ind: individual to be processed
    \param N: the neuron to be deleted
  */
  void delNeuron(Individual& Ind, unsigned N);

  /*__ Ein verstecktes Neuron wird eingefuegt. Dieser erfolgt an der
     angegebenen Stelle, wobei die moeglichen freien Stellen von 0
     beginnend gezaehlt werden. */
  //! Insert a hidden neuron according to a given number.
  /*! All hidden neurons are enumerated so that an arbitrary hidden neuron
    can be determined by an integer number.
    \param Ind: individual to be processed
    \param N: the neuron to be inserted
  */
  void insNeuron(Individual& Ind, unsigned N);

  /*__ Funktion realisiert ein Gluecksrad und liefert den Index des
     erdrehten Feldes (beginnend mit 0) zurueck. Die relativen
     Groessen und die Anzahl der Felder wird durch das Array
     bestimmt. Die hier angegebenen Zahlen muessen nicht auf 1
     normiert sein. */
  //! randomly select a position in a probability matrix and return this position
  /*!
    \param prob: probability matrix
  */

  unsigned wheelOfFortune(Array<double>& prob);
};

#endif

#ifndef FFMANIP_H
#define FFMANIP_H

#include <Array/Array.h>
#include <Array/ArrayIo.h>
#include <Rng/GlobalRng.h>

//#pragma interface
using namespace std;
//! This class paves the way for mutating the struture of neural networks
/*! Operations include counting inputs and outputs of a given neuron,
 checking which neurons or connections could be deleted, and where
neurons and connections could be added.
*/
class FFManip{
 protected:
  //! initialize the connections randomly of a block of non-layered feedforward neural network
  /*!
    \param nIn: number of input nodes
    \param nHid: number of hidden nodes
    \param nOut: number of output nodes
    \param prob: probability for a connection
  */
  Array<int> initFFBulk(const unsigned nIn, const unsigned nHid, const unsigned nOut, const double prob); 

  //! initialize the connections randomly of a block of layered feedforward neural network
  /*!
    \param nIn: number of input nodes
    \param nHid: number of hidden nodes
    \param nOut: number of output nodes
    \param prob1: probability for a connection between input and hidden nodes
    \param prob2: probability for a connection between hidden and output nodes
    \param prob3: probability for a connection between bias nodes and other nodes
  */
  Array<int> initFFOL(const unsigned nIn, const unsigned nHid, 
		      const unsigned nOut, const double prob1, 
		      const double prob2, const double prob3);
//! Check if a given neuron structure is valid, otherwise repair it
/*! 
  \param c: a connection matrix to be checked
  \param nIn: number of input nodes
  \param nHid: number of hidden nodes
  \param nOut: number of output nodes
*/
  void repairBulk(Array<int>& c, const unsigned nIn, const unsigned nHid, 
  const unsigned nOut);
  //! Check if a given connection matrix for a layered NN structure is valid
  /*! 
    \param c: a connection matrix to be checked
    \param nIn: number of input nodes
    \param nHid: number of hidden nodes
    \param nOut: number of output nodes
  */ 
  void repairOL(Array<int> & c, const unsigned nIn, const unsigned nHid, 
  const unsigned nOut);
  //! Add a connection from a given neuron to a randomly selected neuron
  /*! The weight for the newly added connection is generated randomly within
a given interval
    \param c: connection matrix
    \param w: weight matrix
    \param pos: the neuron from which a new connection should be added
    \param low: lower bound of the interval for generating a random weight
    \param high: upper bound of the interval for generating a random weight
  */
  void connectNew(Array<int> & c, Array<double>& w, unsigned pos, double low, double high);
  //! Update connections according to a given matrix to a current one and generate a random weight for that connection
  /*! 
    \param cOld: the old connection matrix
    \param c: the new connection matrix 
    \param w: the weight matrix
    \param low: lower bound of the interval for generating a random weight
    \param high: upper bound of the interval for generating a random weight
  */
  void pullThrough(Array<int>& cOld, Array<int>& c, Array<double>& w, unsigned, double low, double high);
 protected:
  /*__ Berechnet das Maximum der Argumente */
  //! Return \c max(a,b)

  unsigned max(unsigned a, unsigned b);

  /*__ Berechnet das Minimum der Argumente */
  //! Return \c min(a,b)
  unsigned min(unsigned a, unsigned b);

  /*__ Diese Methode zaehlt die Anzahl der Eingaenge des Neurons an
       der uebergebenen Position in der Verbindungsmatrix. Hierbei
       wird der Bias nicht als Eingang mitgezaehlt. */
  //! Count the number of inputs to a given neuron
  /*!
    \param c: connection matrix
    \param N: the neuron whose inputs are to be counted
  */
  unsigned nInput(Array<int>& c, const unsigned N);

  /*__ Diese Methode zaehlt die Zahl der Ausgaenge eines gegebenen
       Neurons. */
  //! Count the number of outputs to a given neuron
  /*!
    \param c: connection matrix
    \param N: the neuron whose outputs are to be counted
  */
  unsigned nOutput(Array<int>& c, const unsigned N);

  /*__ Diese Methode erstellt eine Liste aller moeglicherweise
       loeschbarer Verbindung unter der Randbedingung, dass das so
       erzeugte Netz gueltig bleibt. Somit werden nur solche
       Verbindungen als loeschbar vorgeschlagen, deren Anfangs- und
       Endneurone immernoch nach dem Loeschen mindestens einen Ausgang
       bzw. Eingang behalten. Unter der Annahme, dass nur gueltige
       Netze betrachtet werden, bleiben die Netze auch gueltig. Auserdem werden 
       alle bestehenden Bias-Verbindungen zum Loeschen vorgeschlagen. */
  //! Check which connections can be deleted
  /*! A connection can be removed if after deletion, all neurons have at least
    one input and one output
    \param c: connection matrix
    \param nIn: number of input nodes
    \param nHid: number of hidden nodes
  */
 Array<unsigned> dispensableFFConnection(Array<int>& c, const unsigned nIn, const unsigned nHid);

  /*__ Diese Methode erzeugt eine Liste aller freien Stellen im
       Netzwerk, welche durch eine zusaetzliche Verbindung gefuellt
       werden koennten. Neben der Verbindungsmatrix ist die Anzahl der
       Eingaenge und versteckten Neurone zu uebergeben. */
  //! Check where a connection can be added
  /*! 
    \param c: connection matrix
    \param nIn: number of input nodes
    \param nHid: number of hidden nodes
  */
  Array<unsigned> feasibleFFConnection(Array<int>& c, const unsigned nIn, const unsigned nHid);
};
#endif

Regression and Classification Models Toolbox.
To achieve most flexibility, the ReClaM package is build like a
construction kit, i.e. a network is divided into three major
components: The structure and the behavior of the network, named
here as "model", the algorithm used for the optimization of the
free parameters (weights) of the network and the error measures
that are used by the optimization algorithms to minimize the
network error. You can choose predefined types of each of the
components to build your desired network or add new types by your
own. To manage the communication between the different components,
the class Model Interface is used. 
This library was originally developped at the Institut fuer
Neuroinformatik, Ruhr-Universitaet Bochum, and released under GPL.

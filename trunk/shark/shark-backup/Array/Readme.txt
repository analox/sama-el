The class array serves as a basic data structure for building arrays of arbitrary types
and dimensions.
It is implemented via a C++ template class array< type > and it is possible to dynamically
resize or change the number of dimensions.
Furthermore, parameter passing by value and assignment works properly
(i.e., the value is passed or assigned and not a pointer to the value) and the subscript
operator [ ] may perform a range check at run-time.

//===========================================================================
/*!
 *  \file ArrayCheck.h
 *
 *  \brief Used for exception handling or simulation of exception
 *         handling, when older compilers are used.
 *
 *  The preprocessor functions and the class defined in this file were created
 *  at a time where exceptions were part of the C++ language
 *  only theoretically. <br>
 *  When using older compilers it is possible to simulate exception
 *  handling. The thrown exception is then written to the standard
 *  error stream and the program is interrupted. <br>
 *  For this simulation the flag #HANDLE_EXCEPTIONS below must be
 *  set undefined.<br>
 *  Newer compilers now allow the full usage of exceptions, so
 *  when setting flag #HANDLE_EXCEPTIONS below "defined", you will not 
 *  receive an error, but the running program will exit without any 
 *  specific error message. 
 *  It is then your task to handle any exceptions by using the
 *  C++ try-catch-construct.  
 *
 *  \author  M. Kreutz
 *  \date    1995
 *
 *  \par Copyright (c) 1995, 1999:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *
 *  \par Project:
 *      Array
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: ArrayCheck.h,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:21 $
 *
 *  \par Changes:
 *      $Log: ArrayCheck.h,v $
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.6  2002/05/16 13:18:34  rudi
 *      doxygen commands added/modified.
 *
 *
 *  This file is part of Array. This library is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 2, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *
 */
//===========================================================================

#ifndef __ARRAYCHECK_H
#define __ARRAYCHECK_H

#include <iostream>

#ifndef NDEBUG
#  ifndef CHECK_ASSERTIONS
#    define CHECK_ASSERTIONS
#  endif
#  ifndef CHECK_RANGES
#    define CHECK_RANGES
#  endif
#  ifndef CHECK_SIZES
#    define CHECK_SIZES
#  endif
#  ifndef CHECK_TYPES
#    define CHECK_TYPES
#  endif
#  ifndef CHECK_IO
#    define CHECK_IO
#  endif
     //! If this flag is defined, all exceptions must be handled
     //! by the try-catch-concept of C++ (newer compilers). <br>
     //! If this flag is undefined, then throwing exceptions
     //! is simulated, by writing the exception type
     //! to the standard error stream.
#    define HANDLE_EXCEPTIONS
#endif

#ifndef __EXCEPTION
#define __EXCEPTION

//===========================================================================
/*!
 *  \brief Exception class, that can be used to handle exceptions.
 *         
 *  When an exception is thrown, this class offers methods to store
 *  the type of occurred exception and to return it, when the
 *  exception is handled.
 *
 *  \par Example
 *  \code
 *  ...
 *  try 
 *  {
 *      // Code, where check exceptions can occurr:
 *      ...
 *  }
 *  // Handle the exceptions:
 *  catch ( check_exception &check ) 
 *  {
 *      cout << "Error: Check exception thrown!" << endl;
 *      // Extract check exception type:
 *      cout << "Exception type: " << check.what() << endl;
 *      exit( 8 );
 *  }
 *  \endcode
 *  
 *  \author  M. Kreutz
 *  \date    1995
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class check_exception {
private:
  //! The type of the last called check exception. 
  const char* desc;
public:
 //========================================================================
 /*!
  *  \brief Will store the check exception type "what_arg".
  *
  *  \em what_arg is stored in #desc, so it can be queried by method
  *  #what later.
  *
  *  \param what_arg the check exception type
  *  \return none
  *
  *  \author  M. Kreutz
  *  \date    1995
  *
  *  \par Changes
  *      none
  *
  *  \par Status
  *      stable
  *
  */
  check_exception( const char* what_arg ) : desc( what_arg ) { }

 //========================================================================
 /*!
  *  \brief Will return the stored check exception type.
  *
  *  Will return what was previously stored in #desc by method
  *  #check_exception.
  *
  *  \return the stored check exception type
  *
  *  \author  M. Kreutz
  *  \date    1995
  *
  *  \par Changes
  *      none
  *
  *  \par Status
  *      stable
  *
  */
  const char* what() const { return desc; }
};
#endif

#ifndef RANGE_CHECK
#  ifdef CHECK_RANGES
   //========================================================================
   /*!
    *  \brief Will throw an "range check error" check exception when condition
    *         "cond" is violated.
    *
    *  Normally \em cond is a comparison like \f$i < j\f$, so if
    *  \f$i \geq j\f$ then the exception is thrown.
    *
    *  \param cond the condition that must be fulfilled to prevent a
    *              check exception.
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    1995
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
#    define RANGE_CHECK( cond ) if( ! ( cond ) ) THROW( check_exception( "range check error" ) );
#  else
#    define RANGE_CHECK( cond )
#  endif
#endif

#ifndef SIZE_CHECK
#  ifdef CHECK_SIZES
   //========================================================================
   /*!
    *  \brief Will throw an "size mismatch" check exception when condition
    *         "cond" is violated.
    *
    *  Normally \em cond is a comparison like \f$i < j\f$, so if
    *  \f$i \geq j\f$ then the exception is thrown.
    *
    *  \param cond the condition that must be fulfilled to prevent a
    *              check exception.
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    1995
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
#    define SIZE_CHECK( cond ) if( ! ( cond ) ) THROW( check_exception( "size mismatch" ) );
#  else
#    define SIZE_CHECK( cond )
#  endif
#endif

#ifndef TYPE_CHECK
#  ifdef CHECK_TYPES
   //========================================================================
   /*!
    *  \brief Will throw an "type mismatch" check exception when condition
    *         "cond" is violated.
    *
    *  Normally \em cond is a comparison like \f$i < j\f$, so if
    *  \f$i \geq j\f$ then the exception is thrown.
    *
    *  \param cond the condition that must be fulfilled to prevent a
    *              check exception.
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    1995
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
#    define TYPE_CHECK( cond ) if( ! ( cond ) ) THROW( check_exception( "type mismatch" ) );
#  else
#    define TYPE_CHECK( cond )
#  endif
#endif

#ifndef IO_CHECK
#  ifdef CHECK_IO
   //========================================================================
   /*!
    *  \brief Will throw an "I/O error" check exception when condition
    *         "cond" is violated.
    *
    *  Normally \em cond is a comparison like \f$i < j\f$, so if
    *  \f$i \geq j\f$ then the exception is thrown.
    *
    *  \param cond the condition that must be fulfilled to prevent a
    *              check exception.
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    1995
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
#    define IO_CHECK( cond ) if( ! ( cond ) ) THROW( check_exception( "I/O error" ) );
#  else
#    define IO_CHECK( cond )
#  endif
#endif

#include <csignal>

#ifndef THROW
#  ifdef HANDLE_EXCEPTIONS
   //========================================================================
   /*!
    *  \brief Will throw a check exception of type "except".
    *
    *  If you use older C++ compilers without full support of
    *  exceptions and you have set the flag #HANDLE_EXCEPTIONS,
    *  then program will be interrupted and the check exception
    *  type \em except will be written to the standard error stream. <br>
    *  If you use newer compilers with full exception support, then
    *  the exception \em except will be thrown directly and you have
    *  to handle it and extract its type by your own, using the
    *  try-catch-construct. 
    *
    *  \param except type of the thrown check exception
    *
    *  \author  M. Kreutz
    *  \date    1995
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
#    define THROW( except ) { throw except; }
#  else
#    include <iostream>
#    ifdef __GNUC__
#        define THROW( except ) { std::cerr << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ << ": " << except.what( ) << std::endl; raise( SIGQUIT ); }
#    else
#        define THROW( except ) { std::cerr << __FILE__ << ":" << __LINE__ << ": " << except.what( ) << std::endl; exit( 1 ); }
#    endif
#  endif
#endif


//========================================================================
/*!
 *  \brief Will throw an "undefined operator" check exception when used.
 *
 *  \return none
 *
 *  \author  M. Kreutz
 *  \date    1995
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
#define UNDEFINED { THROW( check_exception( "undefined operator" ) ); }

#endif /* !__ARRAYCHECK_H */








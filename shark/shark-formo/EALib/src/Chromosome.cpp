/*! Chromosome.cpp
* ======================================================================
*
*  File    :  Chromosome.cpp
*  Created :  1998-08-11
*
*  Copyright (c) 1995-2000 Martin Kreutz
*
*  Institut fuer Neuroinformatik
*  Ruhr-Universitaet Bochum
*  44780 Bochum, Germany<BR>
*  Phone: +49-234-32-25558<BR>
*  Fax:   +49-234-32-14209<BR>
*  eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
*  www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
*  
*    \par Project:
*        EALib
*  
*    \par Language and Compiler:
*        C++, egcs (Linux)
*  
*    \par File and Revision:
*        $RCSfile: Chromosome.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: Chromosome.cpp,v $
*        Revision 2.3  2004/06/17 23:08:59  saviapbe
*        An error in the INT's URL was corrected.
*
*        Revision 2.2  2004/06/17 22:53:19  saviapbe
*        The standard file header was for doxygen adapted.
*
*        Revision 2.1  2004/05/20 12:38:05  shark-admin
*        Chromosome::pvm_pkchrom() and Chromosome::pvm_upkchrom() default dummy methods added
*  
*        Revision 2.0  2003/11/28 16:23:09  shark-admin
*        Revision tag reset to revision tag 2.x
*  
*        Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
*        INI Administration
*  <BR>
*  
*  
*  Last update: 20.11.2002 by Marc Toussaint and Stefan Wiegand (INI)
*
* ----------------------------------------------------------------------
*
*  This file is part of the EALib. This library is free software;
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
* ======================================================================
*/

#include "Array/ArrayCheck.h"
#include "EALib/Chromosome.h"

bool Chromosome::operator == ( const Chromosome& c ) const
{
    UNDEFINED
    return false;
}

bool Chromosome::operator < ( const Chromosome& c ) const
{
    UNDEFINED
    return false;
}

// automatically generated
/*
bool Chromosome::operator != ( const Chromosome& c ) const
{
    return !( *this == c );
}

bool Chromosome::operator > ( const Chromosome& c ) const
{
    return c < *this;
}

bool Chromosome::operator <= ( const Chromosome& c ) const
{
    return !( c < *this );
}

bool Chromosome::operator >= ( const Chromosome& c ) const
{
    return !( *this < c );
}
*/

//===========================================================================

//
// added by Marc Toussaint and Stefan Wiegand at 20.11.2002
//

/*! inteface methods for more externally from the EALib defined chromosomes */
void Chromosome::init(){}
void Chromosome::init(const char* filename){}
void Chromosome::mutate(){}
void Chromosome::registerIndividual(const Individual& i,uint you){}
void Chromosome::appendToIndividual(Individual& i){}

/*! PVM routines */
int  Chromosome::pvm_pkchrom()
  { 
    std::cerr << "EALib/Chromosome.cpp: default dummy routine for pvm_pkchrom() implemented." << std::endl;
    return -1 ;
  }
int  Chromosome::pvm_upkchrom()
  {
    std::cerr << "EALib/Chromosome.cpp: default dummy routine for pvm_upkchrom() implemented." << std::endl;
    return -1;
  }

//===========================================================================
//

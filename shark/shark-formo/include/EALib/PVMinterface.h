//===========================================================================
/*!
 *  \file PVMinterface.h
 *
 *  \brief This file serves as a dummy interface for the use of PVM with the EALib.
 *
 *
 *  \author  S. Wiegand
 *  \date    2003-01-01
 *
 *  \par Copyright (c) 1999-2003:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 
 *  \par Project:
 *      EALib
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: PVMinterface.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: PVMinterface.h,v $
 *      Revision 2.3  2004/05/20 12:36:30  shark-admin
 *       #include <stdio.h> added to prevent a compilation error caused by pvm3.h, pvm_pkbyte(...) and pvm_pkbyte(...) added in order to allow the work with variables of type char.
 *
 *      Revision 2.2  2004/03/05 16:32:50  shark-admin
 *      include file iostream added
 *
 *      Revision 2.1  2003/12/02 17:44:19  shark-admin
 *      Bug with autoconf's PVM-check resolved
 *      An additional autoconf-test for Lesstif-check added
 *      Bug with PVM-Interface resolved
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *<BR>
 *
 *
 *  This file is part of EALib. This library is free software;
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

#ifndef PVMINTERFACE_H
#define PVMINTERFACE_H

#define PvmDataDefault  0 

/*! This file serves as a dummy interface for the use of PVM with the EALib. It
    reimplements basic PVM methods with default return error values due to calls 
    implemented in Individual.h. 
*/

typedef unsigned uint;

#ifdef PVM_EXISTS
#include <stdio.h>
#include <pvm3.h>
#else
#include <iostream>
#ifdef __cplusplus
extern "C" {
#endif

int pvm_joingroup(char*          );
int pvm_mytid    (               ); 
int pvm_gettid   (char*  ,int    );

int pvm_initsend (int            );
int pvm_send     (int    ,int    );
int pvm_recv     (int    ,int    );

int pvm_pkdouble (double*,int,int);
int pvm_pkuint   (uint*  ,int,int);
int pvm_pkint    (int*   ,int,int);
int pvm_pkbyte    (char*   ,int,int);

int pvm_upkdouble(double*,int,int);
int pvm_upkuint  (uint*  ,int,int);
int pvm_upkint   (int*   ,int,int);
int pvm_upkbyte   (char*   ,int,int);

int pvm_probe    (int    ,int    );
int pvm_barrier  (char*  ,int    );
int pvm_exit     (               );

int pvm_getinst  (char*  ,int    );
int pvm_bufinfo  (int, int*, int*, int*);

#ifdef __cplusplus
}
#endif

#endif

#endif

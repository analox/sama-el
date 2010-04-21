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
 *      <BR>

 *  \par Project:
 *      EALib
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: PVMinterface.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: PVMinterface.cpp,v $
 *      Revision 2.3  2004/05/20 12:39:13  shark-admin
 *      pvm_pkbyte(...) and pvm_upkbyte(...) default dummy methods defined
 *
 *      Revision 2.2  2004/03/05 16:31:55  shark-admin
 *      some warnings added
 *
 *      Revision 2.1  2003/12/02 17:44:19  shark-admin
 *      Bug with autoconf's PVM-check resolved
 *      An additional autoconf-test for Lesstif-check added
 *      Bug with PVM-Interface resolved
 *
 *      Revision 2.0  2003/11/28 16:23:09  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
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

#include "EALib/PVMinterface.h"

#ifndef PVM_EXISTS

int pvm_joingroup(char* group                          )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_mytid    (                                     )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_gettid   (char* group ,int inum                )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_initsend (int encoding                         )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_send     (int tid     ,int msgtag              )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_recv     (int tid     ,int msgtag              )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_pkdouble (double *dp  ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_pkuint   (unsigned int *up    ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_pkint    (int *ip     ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_pkbyte    (char *xp     ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_upkdouble(double *dp  ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_upkuint  (unsigned int *up    ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_upkint   (int *ip     ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_upkbyte    (char *xp     ,int nitem ,int stride   )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};


int pvm_probe    (int tid     , int msgtag             )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_barrier  (char* group ,int count               )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_exit     (                                     )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_getinst  (char* group , int tid                )
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

int pvm_bufinfo (int bufid, int *bytes, int *mstag, int *tid)
{ 
  std::cerr << "WARNING: Using EALib/PVMinterface" << std::endl; 
  return -1; 
};

#endif

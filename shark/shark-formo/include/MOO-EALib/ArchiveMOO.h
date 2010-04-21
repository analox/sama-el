/* ======================================================================
 *
 *  File    :  ArchiveMOO.h
 *  Created :  2002-01-23
 *
 *  Description : Class ArchiveMOO in EALib for MOO
 *  Author  : Tatsuya Okabe <tatsuya.okabe@honda-ri.de>
 *  Copyright (c) 2002-2004 GPL2 Honda Research Institute Europe GmbH
 * 
 *  \par Maintained by:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>

 *  \par Project:
 *      MOO-EALib
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: ArchiveMOO.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: ArchiveMOO.h,v $
 *      Revision 2.6  2004/12/30 14:15:04  shark-admin
 *      ; added
 *
 *      Revision 2.5  2004/12/30 11:52:48  shark-admin
 *      saveArchiveGPT() added: resulting file is readable for gnuplot
 *
 *      Revision 2.4  2004/11/15 15:36:17  shark-admin
 *      compatibility with vc++ ensured
 *
 *      Revision 2.3  2004/03/03 13:23:49  saviapbe
 *      Copyright information was changed.
 *
 *      Revision 2.2  2004/02/12 11:27:09  saviapbe
 *      Author's contact address, INI Header are modified.
 *
 *      Revision 2.0  2003/11/28 16:23:11  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *<BR>
 
 * ----------------------------------------------------------------------
 *
 *  This file is part of the MOO-EALib. This library is free software;
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
 * ====================================================================== */
 //
 // 	Authors message
 //======================================================================
 /*	Thank you very much for your interest to MOO-EALib.

	Since our company's name was changed on 1st, January, 2003,
	my E-mail address in the source codes were also changed.
	The current E-mail address (6th,Feb.,2004) is as follows:

	tatsuya.okabe@honda-ri.de.

	If you cannot contact me with the above E-mail address,
	you can also use the following E-mail address:

	t_okabe_de@hotmail.com.

	If you have any questions, please don't hesitate to
	ask me. It's my pleasure.

	Best Regards,
	Tatsuya Okabe

	*********************************************************
	Tatsuya Okabe
	Honda Research Institute Europe GmbH
	Carl-Legien-Strasse 30, 63073 Offenbach/Main, Germany
	Tel: +49-69-89011-745
	Fax: +49-69-89011-749
	**********************************************************/


#ifndef __ARCHIVEMOO_H
#define __ARCHIVEMOO_H

#include "MOO-EALib/PopulationMOO.h"
#include "MOO-EALib/IndividualMOO.h"
#include "Array/Array.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "fstream"
#ifdef _WIN32
 #include <limits>
#else
 #include <values.h>
#endif
//#include "values.h"

class ArchiveMOO : private std::vector< IndividualMOO * >
{
 public:
  //*******************************************
  //** Constructor
  //*******************************************
  //***** TO-AM-001
  ArchiveMOO( );
  //***** TO-AM-002
  ArchiveMOO( bool strategy );
  //***** TO-AM-003
  ArchiveMOO( unsigned max );
  //***** TO-AM-004
  ArchiveMOO( unsigned max, bool strategy );

  //*******************************************
  //** Destructor
  //*******************************************
  //***** TO-AM-005
  ~ArchiveMOO( );

  //*******************************************
  //** Internal variables I/O
  //*******************************************
  //***** TO-AM-010
  unsigned getMaxArchive( );
  //***** TO-AM-011
  void     setMaxArchive( unsigned max );
  //***** TO-AM-012
  unsigned getCapacity( );
  //***** TO-AM-013
  unsigned size( );
  //***** TO-AM-014
  bool     getStrategy( );
  //***** TO-AM-015
  void     setStrategy( bool strategy );
  //***** TO-AM-016
  void     minimize( );
  //***** TO-AM-017
  void     maximize( );
  
  //*******************************************
  //** Archive I/O
  //*******************************************  
  //***** TO-AM-050
  void addArchive( IndividualMOO& indmoo );
  //***** TO-AM-051
  IndividualMOO& readArchive( unsigned i );
  //***** TO-AM-052
  void delArchive( unsigned i );
  //***** TO-AM-053
  void delArchive( std::vector< unsigned >& v );
  //***** TO-AM-054
  void delArchive( Array< unsigned >& a );
  //***** TO-AM-055
  void cleanArchive( );
  //***** TO-AM-056
  void delSharingWorst( );
  //***** TO-AM-057
  void delSharingWorst( double div );
  //***** TO-AM-058
  IndividualMOO& readBestArchive( );
  //***** TO-AM-059
  void nonDominatedSolutions( );
  
  //*******************************************
  //** Dominate
  //*******************************************  
  //***** TO-AM-100
  int Dominate( IndividualMOO& im1, IndividualMOO& im2 );
  //***** TO-AM-101
  int Dominate( IndividualMOO& im1 );
  //***** TO-AM-102
  void delDominateArchive( IndividualMOO& im1 );

  
  //*******************************************
  //** Distance on Fitness Space
  //*******************************************  
  //***** TO-AM-150
  double distanceOnFitness( unsigned i1, unsigned i2 );
  //***** TO-AM-151
  double distanceOnFitness( IndividualMOO& im1 );
  //***** TO-AM-152
  Array< double > distanceDataOnFitness( );
  //***** TO-AM-153
  unsigned sharingWorst( );
  //***** TO-AM-154
  unsigned sharingWorst( double div );
  //***** TO-AM-155
  unsigned sharingBest( );
  //***** TO-AM-156
  double minDistanceOnFitness( );
  //***** TO-AM-157
  double minDistanceOnFitness( unsigned i1 );
  

  //*******************************************
  //** PAES
  //*******************************************
  //***** TO-AM-300
  unsigned crowded( IndividualMOO& im, double div );
  


  //*******************************************
  //** For Metrix
  //*******************************************  
  //***** TO-AM-1000
  void saveArchive( char *filename );

  //***** SW-AM-1001
  void saveArchiveGPT( char *filename );



 protected:
  unsigned MaxArchive;
  bool     Strategy;



};

#endif /* !__ARCHIVEMOO_H */

/* ZNminmax.h */
/* ======================================================================
 *
 *  File    :  ZNminmax.h
 *  Created :  2001-08-17
 *
 *  Copyright (c) 1995-2001 Martin Kreutz
 *  EMail - Martin.Kreutz@zn-ag.de
 *
 *  ZN Vision Technologies AG
 *  Universitaetsstrasse 160
 *  44801 Bochum
 *  Germany
 *  
 *  Last update: 2001-08-17
 *
 * ----------------------------------------------------------------------
 *
 *  This file is part of the Shark distribution. Shark is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 2, or (at your option) any later version.
 *
 *  Shark is distributed in the hope that it will be useful,
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


/* 
 *
 * As the new Shark-version (1.0.8) was arranged at the ZN 
 * (Zentrum fuer Neuroinformatik, see above),
 * it was necessary to respect the general directory structure used
 * there. The file "ZNminmax.h" is used for the compilation
 * with the Visual C++ compiler.
 * The original file includes another file with the real definitions,
 * so all was put together to form this dummy file to spread with
 * the "non-ZN-version".
 *
 */ 


#ifndef ZN_MINMAX_H
#define ZN_MINMAX_H

#ifdef	_MSC_VER

#ifndef	NOMINMAX
#define NOMINMAX // prevent Microsoft from defining min and max in windef.h
#endif
#if _MSC_VER < 1300   // because .NET has owned macros with the same names and causes an error
namespace std
{

	
template <class T>
inline const T& min(const T& a, const T& b)
{
	return b < a ? b : a;
}

template <class T>
inline const T& max(const T& a, const T& b)
{
	return a < b ? b : a;
}

template <class T, class Compare>
inline const T& min(const T& a, const T& b, Compare comp)
{
	return compare(b, a) ? b : a;
}

template <class T, class Compare>
inline const T& max(const T& a, const T& b, Compare comp)
{
	return compare(a, b) ? b : a;
}

}
#endif // _MSC_VER 
#else
#ifdef	__GNUC__ 
#		include <algorithm>
#else
#		error "cannot determine the system and/or compiler type"
#endif
#endif

#endif // ZN_MINMAX_H



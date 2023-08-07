/*=============================== extras ===============================*/
/*

   Additional definitions, constants, and macros that are defined
   in case simple implementations change from one platform to the next.


   Name:	extras.h

   Author:	Patricio A. Vela, pvela@ece.gatech.edu

   Created:	08/XX/2004
   Modified:	11/09/2004

*/
/*=============================== extras ===============================*/

#ifndef _EXTRAS_H
#define _EXTRAS_H

#ifndef NULL
#define NULL	0
#endif

#define ISNULL(ptr)	!(ptr)
#define NOTNULL(ptr)	(ptr)

#define STRINGIFY(ARG)	#ARG
#define TOSTRING(ARG)	STRINGIFY(ARG)


#if defined (_WIN64_)
#  define ISNAN(a)	((a) != (a))
#  define ISINF(x)	!(_finite(x))
#else 
#  if defined(WIN32) && defined(_MSC_VER)
#    define ISNAN		_isnan
#    define ISINF(x)	!(_finite(x))
#  else
#    define ISNAN	isnan
#    define ISINF	isinf
#  endif
#endif
 


#endif /* _EXTRAS_H */

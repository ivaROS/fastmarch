/*=============================== imtype ===============================*/
/*

  Name:		imtype.h

  Athor:	Patricio A. Vela, pvela@ece.gatech.edu

  Created:	01/18/2006
  Modified:	01/18/2006


  Header file for an image type implementation.  This is so that different
  C/C++ files of the same project can have a uniform definition of what an
  image is, but which can be modified through compile time options.

*/
/*=============================== imtype ===============================*/
#ifndef __IMTYPE_H
#define __IMTYPE_H


#if defined(_IMFLOAT_)

typedef float imtype;

#ifdef _MEX_
#define mxIMAGE_CLASS mxSINGLE_CLASS
#endif

#elif defined(_IMDOUBLE_)

typedef double imtype;

#ifdef _MEX_
#define mxIMAGE_CLASS mxDOUBLE_CLASS
#endif

#elif defined(_IMINTEGER_)
typedef int imtype;

#ifdef _MEX_
#define mxIMAGE_CLASS mxINT32_CLASS
#endif

#else

typedef unsigned char imtype;

#ifdef _MEX_
#define mxIMAGE_CLASS mxUINT8_CLASS
#endif

#endif

typedef double imbinary;
#ifdef _MEX_
#define mxIMBINARY_CLASS mxDOUBLE_CLASS
#endif

#endif

/*
*/
/*=============================== imtype ===============================*/

/*============================== number.h ==============================*/
/*
   Name:	number.h

   Author:	Patricio A. Vela, pvela@ece.gatech.edu

   Created:	2004/04/05
   Modified:	2008/09/12


   Header file for implementation of a generic number type.  A number
   can be integer based or real number based, so long as the basic 
   arithmetic operations are defined for it.


   History:

   ***** v0.1 *****
   *
   * -created.
   *
   *****
*/
/*============================== number.h ==============================*/

#ifndef _NUMBER_H
#define _NUMBER_H


#if defined(_INTEGER_)

/*--Jimi is going to help me redo this.--*/
typedef int number;
#define ISHIFT		8
#define MAXNUM 		(1<<(2*ISHIFT))
#define TIMES(x,y)	(((x)*(y))>>ISHIFT)
#define DIVIDE(x,y)	(((x)/(y))<<ISHIFT)
#define POW2(x,n)	((x)<<(n))
#define IPOW2(x,n)	((x)>>(n))
#define TIM2(x)		((x)<<1)
#define DIV2(x)  	((x)>>1)
#define INV(x)   	((number)((1<<ISHIFT)/(x)))
#define NCONV(x)	(((number)(x))<<ISHIFT)
#define NUMBER(x)	((number)((x)*(1<<ISHIFT)))
#define ABS(x)		abs((int)(x))

#elif defined(_DOUBLE_)

#include<math.h>
typedef double number;
#ifdef _MEX_
#define mxNUMBER_CLASS	mxDOUBLE_CLASS
#endif

#define TIMES(x,y)	((x)*(y))
#define DIVIDE(x,y)	((x)/(y))
#define POW2(x,n)	((x)*((real)(1<<(n))))
#define IPOW2(x,n)	((x)/((real)(1<<(n))))
#define TIM2(x)		(2*(x))
#define DIV2(x)  	((x)/2.0)
#define INV(x)   	(1.0/(x))
#define NCONV(x)	(x)
#define NUMBER(x)	(x)
#define SQRT		sqrt
#define ABS		fabs
#define EXP		exp
#define LOG		log
#define POW		pow
#define ATAN		atan
#  if defined(_WIN64_) || defined(_WIN32) || defined(_MSC_VER)
 #define ROUND(x) ((x-floor(x))>0.5 ? ceil(x) : floor(x))
#  else
#define ROUND		round
#endif

#else

#include<math.h>

typedef float number;
#ifdef _MEX_
#define mxNUMBER_CLASS	mxSINGLE_CLASS
#endif

#define MAXNUM 		(1e10)
#define TIMES(x,y)	((x)*(y))
#define DIVIDE(x,y)	((x)/(y))
#define POW2(x,n)	((x)*((float)(1<<(n))))
#define IPOW2(x,n)	((x)/((float)(1<<(n))))
#define TIM2(x)		(2*(x))
#define DIV2(x)  	((x)/2.0)
#define INV(x)   	(1.0/(x))
#define NCONV(x)	(x)
#define NUMBER(x)	(x)
#define SQRT		sqrtf
#define ABS		fabsf
#define EXP		expf
#define LOG		logf
#define POW		powf
#define ATAN		atanf

#  if defined(_WIN32) || defined(_MSC_VER)
inline double ROUND(double x) { return (x-floor(x))>0.5 ? ceil(x) : floor(x); }
#  else
#define ROUND		roundf
#endif

#endif


#endif /* _NUMBER_H */
/*============================== number.h ==============================*/

/*=============================== mymex ===============================*/
/*

  Author:       Patricio A. Vela, pvela@ece.gatech.edu

  Created:      03/16/2004
  Modified:     05/18/2004


  This file contains a variety of definitions that are typically
  used by code interfacing with MATLAB via the mex interface.
  It also allows for alternative definitions to be used in case
  the mex-compilation is turned off.  To enable mex-interface
  definitions use the -D_MEX_ switch with compilation.

  The mymex header file in the root source directory should be
  though of as a template hear file for use by c/c++ source code.
  If basic definitions and typedefs are to be changed, it is best
  to copy this file into the working directory first as opposed to
  modifying the globally available mymex (template) header file.

  History:

  ***** v0.1a - 05/18/2004 *****
  *
  * -modified to allow for s-function interface description.
  *  The difference is that although we have MATLAB indexing,
  *  MATLAB memory allocation/destruction commands cannot be 
  *  used.
  *
  *****

  ***** v0.1 - 04/30/2004 *****
  *
  * -created.
  *
  *****

*/
/*=============================== mymex ===============================*/

#ifndef __MYMEX_H
#define __MYMEX_H



/*======== Conditional header elements. ========*/

#if defined(_MEX_)	/* _MEX_ { */

#define MALLOC  mxMalloc
#define CALLOC  mxCalloc
#define REALLOC mxRealloc
#define FREE    mxFree
#define PERSIST		mexMakeMemoryPersistent
#define ARRAYPERSIST	mexMakeArrayPersistent
#define ARRAYFREE	mxDestroyArray


#define PRINTF 	mexPrintf		/* Print function. */
#define DEBUGF 	mexPrintf		/* Debug print function. */ 
#define ERRORF 	mexErrMsgIdAndTxt	/* Error print function. */

typedef double real;
typedef int integer;

#define mxCreateRealMatrix	mxCreateDoubleMatrix
#define mxREAL_CLASS		mxDOUBLE_CLASS

#define mexEXP		exp


#elif defined(_SFUNC_)	/* _MEX_ }{ _SFUNC_ */

#define MALLOC  malloc
#define CALLOC  calloc
#define REALLOC realloc
#define FREE    free

#define PRINTF 		ssPrintf		/* Print function. */
#define DEBUGF		ssPrintf		/* Debug print function. */
#define WARNF		ssPrintf		/* Warning print function. */
#define WARNSS(S,msg)	ssWarning(S,msg)	/* Warning print function. */
#define ERRORSS(msg)	ssSetErrorStatus(SS,msg)
#define ERRORF(s1, s2) printf("%s ; %s.\n", s1, s2)
						/* Error print function. */

typedef real_T  real;
typedef int_T integer;

#define mexEXP		exp

#else 			/* _SFUNC_ }{ */

#define MALLOC  malloc
#define CALLOC  calloc
#define REALLOC realloc
#define FREE    free

#define PRINTF	printf			/* Print function. */
#define DEBUGF	printf			/* Debug print function. */
					/* Error print function. */
#define ERRORF(s1, s2) printf("%s ; %s.\n", s1, s2)

typedef float  real;
typedef int integer;

#define mxREAL_CLASS	mxSINGLE_CLASS
#define mexEXP		expf

#endif 			/* } */



/*  In order to use the _MATRIX_ defines below, it is required that 
    NROWS and NCOLS be defined.  This flexibility allows for the macros
    to be used as explicitly needed by the invoking code-space and 
    context.
*/

#ifdef _MATRIX_ /* { */

#if defined(_MEX_) 

#define DX              NROWS
#define DY              1
#define XSUB(ind)       ((int)((ind)/NROWS))
#define YSUB(ind)       ((ind)%NROWS)
#define SUB2IND(x,y)    ((x)*NROWS+(y))
#define MSUB2IND(x,y,nrows,ncols)    ((x)*(nrows)+(y))

#else /* _MEX_ || _SFUNC_  */

#define DX              1
#define DY              NCOLS
#define XSUB(ind)       ((ind)%NCOLS)
#define YSUB(ind)       ((int)((ind)/NCOLS))
#define SUB2IND(x,y)    ((x)+(y)*NCOLS)
#define MSUB2IND(x,y,nrows,ncols)    ((x)+(y)*(ncols))

#endif /* _MEX_ || _SFUNC_  */

#endif /* _MATRIX_ } */


/*========= Permanent header elements. =========*/

typedef unsigned char uint8;

#ifndef _MASKTYPE_
#define _MASKTYPE_ double
#endif
typedef _MASKTYPE_ masktype;


#endif /* __MYMEX_H */

/*=============================== mymex ===============================*/

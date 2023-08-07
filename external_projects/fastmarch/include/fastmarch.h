/*============================== fastmarch =============================*/
/*
   Name:	fastmarch.h

   Author:	Patricio A. Vela, pvela@ece.gatech.edu

   Created:	2004/03/11
   Modified:	2009/03/21

   fastmarch.h is the header file for the fast march implementation.
   

   History:

   ***** v0.1 - 2009/03/21 *****
   *
   * -forked and isolated fast march code from version 2.7 of levelset.
   *
   *****

   Copyright(c) 2005.
*/
/*============================== fastmarch =============================*/

#ifndef __FASTMARCH_H
#define __FASTMARCH_H

/*------------------ Includes ------------------*/
#include<mymex.h>
#include<extras.h>
#include<number.h>
#include<imtype.h>
#include<heap.h>



/*----------------- Definitions ----------------*/
#ifdef __FASTMARCH_CC /* { */

#define ROOT2OVER2	(NUMBER(SQRT(2.0)/2.0))
#define ROOT2		(NUMBER(SQRT(2.0)))

#define MAXNUM		99999
#define INFDIST		90000		/* Inf. distance for level set. */
#define INFTIME		MAXNUM		/* Inf. time for fast marching. */
#define INFSPEED	MAXNUM		/* Inf. speed for ext. velocity. */

#define DEF_MATNAME		"fmdata"	/* Default matlab dump variable name. */
#define DEF_CALLBACK	"marchBack"	/* Default callback name. */

#define maxsq(x,y) ( (x)>(y) ? TIMES(x,x) : TIMES(y,y) )
#define minsq(x,y) ( (x)<(y) ? TIMES(x,x) : TIMES(y,y) )
#define MIN(x,y)   ( (x)<(y) ? (x) : (y) )
#define MAX(x,y)   ( (x)>(y) ? (x) : (y) )

#endif /* __FASTMARCH_CC } */

#define NAMELEN	128		/* Max length of matlab dump variable name.  */

/*--------------- Type Definitions -------------*/
typedef char int8;
#ifndef _DATATYPE_
#define _DATATYPE_ double
#define mxDATA_CLASS mxDOUBLE_CLASS
#endif
#ifndef mxDATA_CLASS
#error  "Need to define mxDATA_CLASS"
#endif
typedef _DATATYPE_ datatype;
typedef int labeltype;
#define mxLABEL_CLASS mxINT32_CLASS

typedef enum{IN=-1,ZERO=0,OUT=1} signtype;	/* Sign values for IN/OUT. */
typedef enum{KNOWN, TRIAL, FAR} fmstat; 	/* Fast marching values. */
typedef enum{ALIVE,QUEUED,SET,DONE}  qstat;	/* Queue status. */
typedef number exptype;
typedef unsigned char shocktype;
#define mxSHOCK_CLASS mxUINT8_CLASS

inline number square(number dV) { return dV*dV; }
inline number max(number dX, number dY) { return (dX>dY?dX:dY); }
inline number min(number dX, number dY) { return (dX<dY?dX:dY); }
inline number pos(number dV) { return (dV>0?dV:0); }
inline number neg(number dV) { return (dV<0?dV:0); }

class fastmarch
{

public:
  fastmarch(void);
  ~fastmarch(void);

  void init( int iM, int iN);
#ifdef _MEX_
  void init(int iM, int iN, char *matname);
  void setName(char *matname);
#endif

  bool checkDimensions(int sM, int sN);

  void setSeedPoint(int x, int y);
  void setSeedPoint(number x, number y);
  void setSeedPoints(int *seeds, int slen);
  void setSeedPoints(int *xpts, int *ypts, int zlen);
  void setSeedRegions(masktype* bwM, int mM, int mN);
  void setKnownRegions(masktype* bwM, int mM, int mN);

  void setDistance( number *sdist );
  void identifyFrontRegion(void);
  void identifyKnownBoundary(void);
  void forceRedistance(void);
  void compdistCostField( number *cost );
  void compdistSpeedField( number *cost );


  number* getDistance(void);
  void getDistance(number *rpsi);

  void setLabelPoints(int *dlabels, int dlen);
  void setLabelPoints(int *xpts,int *ypts, int len);
  void setLabelRegions( labeltype *nlabel, int mM, int mN );
  labeltype* getLabels(void);
  void getLabels(labeltype *rlabels);
  void compLabels(void);
  void compLabelsCostField( number *cost );
  void compLabelsSpeedField( number *cost );
  void compLabelsByCost( number *cost , int len);
  void compLabelsBySpeed( number *cost, int len );

  void compRegionsNMS(number *values, number *thresholds);

  void setData(datatype *data, int iM, int iN, int iP);
  void getData(datatype *data);
  void getDataDimensions(int &iM, int &iN, int &iP);
  bool checkDataDimensions(int sM, int sN, int sP);
  void extendData(void);

  void setPhi(number *nphi, int pM, int pN);
  void addtoPhi(number *aphi);
  void scalePhi(number *sphi);
  number* getPhi(void);

#ifdef _MEX_
  void setParms(mxArray *lsparms);
  void mexPutData(void);
  void mexPersist(void);
#endif

  shocktype* getShockMap(void);
  void getShockMap(shocktype *rshock);

  void free(void);

protected:

  int compZeroContour(void);
  int compZeroContour(number close);
  void compZeroRegion(int* zregion, int& zarea);
  void compZeroRegion(int* zregion, int& zarea, number dist);

  int compKnownBoundary(void);

  void initDist(void);
  void interpDist(void);
  number solveDist(int index, number speed);
  number solveDist(int index, int x, int y, labeltype label, number speed);
  number solveDistL1(int index, number speed);
  number solveDistN8(int index, number speed);
  void distmarch(void);
  void distanceUpdatePDE(void);
  void distanceUpdatePDEHigherOrder(void);
  void redistance(void);
  void distReset(void) { redistance(); }
  void fieldmarch(number *speed);


  void initLabel(void);
  void labelmarch(void);
  void labelmarch(number maxdist);
  void labelfieldmarch(number *speed);
  void labelspeedmarch(number *labelspeed);
  void relabel(void);

  void initNMS(number *values);
  void nmsmarch(number *values, number *thresholds);

  number solveData(int index, number speed);
  void datamarch(void);


  int		M, N, P, 		/* Data scalar and volumetric */
            area, volume;	/*   dimensions.              */

  number	*psi;			/* distance function.         */
  signtype	*sign;			/* sign of distance (in/out). */
  number	*phi;			/* speed/stopping function.   */
  int8		*nogo;			/* forbidden or nogo zones.   */

  labeltype *labels;		/* for propogating labels.    */
  shocktype	*shock;			/* shock map (voronoi).       */
  datatype	*xdata;			/* data proper.               */

  fmstat	*status;		/* status of fast marching.   */
  heap		minheap;		/* heap for fast marching.    */
  int		*seed;			/* the seed data.             */
  int 		seedlen;		/* length of seed list.       */

  #ifdef _MEX_
  bool		isPersistent;		/* Should memory be persisten?    */
  char 		matName[NAMELEN];	/* Name of matlab dump variable.  */
  char		callback[NAMELEN];	/* Name of evolution callback.    */
  mxArray	*matVar;		/* Pointer to matlab variable.    */
  mxArray	*matParms;		/* Pointer to external parm defs. */
  #endif 

};


#endif /* __FASTMARCH_H */
/*============================== fastmarch =============================*/

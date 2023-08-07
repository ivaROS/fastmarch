/*============================== fastmarch =============================*/
/*
   Name:	fastmarch.cc

   Author:	Patricio A. Vela,	pvela@gatech.edu

   Created:	2004/03/11
   Modified:	2009/05/21


   fastmarch is an implementation of fast marching for the Eikonal equation.
   Hopefully the code is sufficiently flexible to allow for a variety of
   implementations with no modifications to the source code (done
   through switches/defines). 


   History:

   ***** v0.1 - 2009/03/21 *****
   *
   * -forked and isolated fast march code from version 2.7 of levelset.
   * -may have retained a few things (like iterative PDE code).
   *
   *****

   Copyright(c) 2005.
*/
/*============================== fastmarch =============================*/

#define __FASTMARCH_CC

/*------------------ includes ------------------*/
#include<math.h>

#if defined(WIN32) && defined(_MSC_VER)
#include<float.h>
#endif

#include<stdlib.h>
#include<stdio.h>
#include<string.h>		/* for memcpy */

#ifdef _MEX_
  #include<mex.h>
#endif

#include<fastmarch.h>


/*------------------- defines ------------------*/

#define DEBUG 1

#ifndef NULL
#  define NULL 0
#endif

#define PI 3.14159


/*================================ code ================================*/


/*----------------------------- fastmarch ----------------------------*/
/*
   Upon instantiation of a fastmarch class object, the constructor 
   initializes all class data to represent nothing.  Parameters are 
   set to default values.

   Inputs: 	N/A.
   Outputs: 	N/A.
*/

fastmarch::fastmarch(void)
{

M = 0;				/* Set to level set with no area. */
N = 0;
P = 0;
area = 0;
volume = 0;

xdata = (datatype *)NULL;		/* Ensure that all pointers are NULL. */
psi  = (number *)NULL;
sign = (signtype *)NULL;
phi  = (number *)NULL;
nogo = (int8*)NULL;

status = (fmstat *)NULL;
seed  = (int *)NULL;
seedlen = 0;

//altPhi = (funcPhi*)NULL;

#ifdef _MEX_
strcpy(matName , DEF_MATNAME);
strcpy(callback, DEF_CALLBACK);
matVar = (mxArray *)NULL;
matParms = (mxArray *)NULL;
#endif

} /*-------------------------- fastmarch --------------------------*/


/*------------------------------- init -------------------------------*/
/*
   Initializes the fastmarch object.  Because this fastmarch 
   implementation is specifically for working on array data, the
   array dimensions are inputs.   Whatever can be initialized with 
   this information is done so (usually that means pre-allocation of 
   working memory).

   Inputs:
     M, N	- array (planar) dimensions.

   Outputs:	an initialized fastmarch object.

*/

#define LDEBUG 0

void fastmarch::init(int iM, int iN) 
{
register int i;

M = iM;			/* Store array dimensions. */
N = iN;
area = iM*iN;

#ifdef _DEBUG_
if (LDEBUG || (DEBUG > 2))
 {
  DEBUGF("Starting fastmarch init.\n");
  DEBUGF("\tImage characteristics: %d, %d -> %d \n", M, N, area);
 }
#endif

/*-- Allocate memory for fast marching structures. --*/
psi  = (number *)   MALLOC(area*sizeof(number));
sign = (signtype *) MALLOC(area*sizeof(signtype));

labels = (labeltype *)NULL;			// Don't use unless needed.
shock = (shocktype *)MALLOC(area*sizeof(shocktype));

status = (fmstat *) MALLOC(area*sizeof(fmstat));
nogo   = (int8 *) MALLOC(area*sizeof(int8));

if ( ISNULL(psi) || ISNULL(sign) || ISNULL(status) 
                 || ISNULL(shock) || ISNULL(nogo) )
 {
  ERRORF("fastmarch::init -- allocProb",
         "Can't allocate memory for fastmarch members.\n");
  return;
 }

#ifdef _DEBUG_
if (LDEBUG & 2)
  DEBUGF("\tDone allocating basic space.\n");
#endif

minheap.init(area);

#ifdef _DEBUG_
if (LDEBUG & 2)
  DEBUGF("\tDone allocating heap space.\n");
#endif

for(i=area-1;i>=0;i--)		/* Initialize ... */
 {
  psi[i]  = INFDIST;		/* Distance function. */
  sign[i] = OUT;		/* Sign function.     */

  status[i] = FAR;		/* Fast march status. */
  nogo[i] = 0;
 }


#ifdef _DEBUG_
if (LDEBUG)
  DEBUGF("Terminating fastmarch init.\n");
#endif

}  

#undef LDEBUG

   /*---------------------------- init ----------------------------*/

#ifdef _MEX_

/*------------------------------- init -------------------------------*/
/*
   Similar functionality as the above init functions, however this
   version is only declared/defined when _MEX_ flag exists.  The 
   function is specific to a MEX implementation.

   Inputs:
     M, N	- array (planar) dimensions.
     matname	- name to set the Matlab output variable to (if used).

   Outputs:	an initialized fastmarch object.

*/

void fastmarch::init(int iM, int iN, char *matname)
{
init(iM, iN);
setName(matname);
}  


   /*---------------------------- init ----------------------------*/


/*------------------------------ setName -----------------------------*/
/*
   Define a new name for the matlab output variable.
*/

void fastmarch::setName(char *newname)
{

strcpy(matName,newname);

}  /*--------------------------- setName --------------------------*/

#endif

/*-------------------------- checkDimensions -------------------------*/
/*
    Returns true if the dimensions do not match, otherwise things are OK.
*/

bool fastmarch::checkDimensions(int sM, int sN)
{
return !( (sM == M)  && (sN == N));
} /*----------------------- checkDimensions ----------------------*/


/*--------------------------- setSeedPoint ---------------------------*/
/*
   Defines the initial (point) front to march out from.  

   Inputs:
     seeds		- indices of seed points. 
     xpt, ypt		- coordinate location of seed.

   Outputs: 		N/A.
   Requirements: 	- seed location is valid.
   Modifies:
     -the zero contour is defined according to the inputs.

   Memory Usage:
     -may allocate space for two integers if input is invalid.
     -allocates sufficient space to store the seed data.

     -for version taking xpts and ypts temporary space is allocated 
      to generate inputs that linear indices version of the setSeeds
      function accepts.

*/

#define NROWS M
#define NCOLS N

//#define LDEBUG 2

void fastmarch::setSeedPoint(int x, int y)
{
int i, j;

if (seedlen > 0)	/* if zero contour exists, wipe it. */
 {
  seedlen = 0;
  ::FREE(seed);
 }

seedlen = 1;
seed = (int *) MALLOC(seedlen*sizeof(int));
if (ISNULL(seed))
 {
  seedlen = 0;
  ERRORF("mexFunction:mallocFail",
         "Unable to allocate memory for Zero Contour.");
  return;
 }

seed[0] = SUB2IND(x,y);		/* set the seed index,      */
j = seed[0];			/* get index value,         */
psi[j]  = 0;			/* set distance to zero,    */
sign[j] = ZERO;			/* set sign to zero, and    */
status[j] = KNOWN;		/* fast march status known. */
}

//NOTE:  May have link issues if number is an int given above function.

void fastmarch::setSeedPoint(number x, number y)
{
int i, j, rx, ry;
int sind;
float t;

if (seedlen > 0)	/* if zero contour exists, wipe it. */
 {
  seedlen = 0;
  ::FREE(seed);
 }

seedlen = 4;
seed = (int *) MALLOC(seedlen*sizeof(int));
if (ISNULL(seed))
 {
  seedlen = 0;
  ERRORF("mexFunction:mallocFail",
         "Unable to allocate memory for Zero Contour.");
  return;
 }

rx = (int)x;			/* round down to int.  */
ry = (int)y;			/* this is top-left corner of 2x2. */

seed[0] = SUB2IND(rx,ry);	/* set the seed index. */
psi[seed[0]] = SQRT((x-rx)*(x-rx) + (y-ry)*(y-ry));
sind = 1;

if ( (rx < (NCOLS-1)) && (ry < (NROWS-1)) )
 {
  seed[sind] = seed[0] + DX;
  if (nogo[seed[sind]] == 0)
   {
    psi[seed[sind]] = SQRT((x-rx-1)*(x-rx-1) + (y-ry)*(y-ry));
    sind++;
   }

  seed[sind] = seed[0] + DX + DY;
  if (nogo[seed[sind]] == 0)
   {
    psi[seed[sind]] = SQRT((x-rx-1)*(x-rx-1) + (y-ry-1)*(y-ry-1));
    sind++;
   }

  seed[sind] = seed[0] + DY;
  if (nogo[seed[sind]] == 0)
   {
    psi[seed[sind]] = SQRT((x-rx)*(x-rx) + (y-ry-1)*(y-ry-1));
    sind++;
   }
  seedlen = sind - 1;
 }
else if (rx < (NCOLS-1))
 {
  seed[sind] = seed[0] + DX;
  if (nogo[seed[sind]] == 0)
   {
    psi[seed[sind]] = SQRT((x-rx-1)*(x-rx-1) + (y-ry)*(y-ry));
	sind++;
   }
  seedlen = sind-1;
 }
else if (ry < (NROWS-1))
 {
  seed[1] = seed[0] + DY;
  if (nogo[seed[sind]] == 0)
   {
    psi[seed[1]] = SQRT((x-rx)*(x-rx) + (y-ry-1)*(y-ry-1));
	sind++;
   }
  seedlen = sind-1;
 }

for (i=seedlen-1;i>=0;i--)	/* for each element of zero contour ... */
 {
  j = seed[i];				/* get index value,         */
  sign[j] = OUT;			/* set sign to out, and    */
  status[j] = KNOWN;		/* fast march status known. */
 }
}

#ifdef LDEBUG
#undef LDEBUG
#endif

#undef NROWS
#undef NCOLS

   /*------------------------ setSeedPoint ------------------------*/

/*--------------------------- setSeedPoints --------------------------*/
/*
   Defines the initial front to march out from.  Typically invoked with 
   a collection of points.  Can be invoked by specifying linear index
   locations or (x,y) image coordinate locations.

   void setSeedPoints( seeds, seed_length);
   void setSeedPoints( x coords, y coords, length);

   Inputs:
     seeds		- indices of seed points. 
     xpts, ypts		- coordinate version of seeds.
     slen		- number of the seed points

   Outputs: 		N/A.
   Requirements: 	N/A.
   Modifies:
     -the zero contour is defined according to the inputs.

   Memory Usage:
     -may allocate space for two integers if input is invalid.
     -allocates sufficient space to store the seed data.

     -for version taking xpts and ypts temporary space is allocated 
      to generate inputs that linear indices version of the setSeeds
      function accepts.

*/

#define NROWS M
#define NCOLS N

//#define LDEBUG 2

void fastmarch::setSeedPoints(int *dseeds, int dslen)
{
int *newptr;
int i, j, x, y;
float t;

newptr = NULL;

if (seedlen > 0)	/* if zero contour exists, wipe it. */
 {
  seedlen = 0;
  ::FREE(seed);
 }

if (dslen <= 0)		/* if passed length is nonsensical ... */
 {
  newptr = (int *)MALLOC(2*sizeof(int));
  if (ISNULL(newptr))
   {
    ERRORF("levelset::setZeroContour::allocFail",
           "Could not allocate temporary memory.");
    return;
   }
  seedlen = 1;		/* set seed point to center of image.  */
  newptr[0] = N>>1;
  newptr[1] = M>>1;
  seed = newptr;
 }
else
 {
  seedlen = dslen;

  seed = (int *) MALLOC(seedlen*sizeof(int));
  if (ISNULL(seed))
   {
    ERRORF("mexFunction:mallocFail",
           "Unable to allocate memory for Zero Contour.");
    return;
   }
 }

memcpy(seed,dseeds,dslen*sizeof(int));


for (i=seedlen-1;i>=0;i--)	/* for each element of zero contour ... */
 {
  x = seed[i];				/* get index value,         */
  psi[x]  = 0;				/* set distance to zero,    */
  sign[x] = ZERO;			/* set sign to zero, and    */
  status[x] = KNOWN;		/* fast march status known. */
 }

} 


void fastmarch::setSeedPoints(int *xpts,int *ypts, int len)
{
int *dseeds, i;

dseeds = (int *)MALLOC(sizeof(int)*len);

if (ISNULL(dseeds))
 {
  ERRORF("levelset::setZeroContour -- mallocFail",
         "Unable to allocate memory for Zero Contour.");
  return;
 }

for (i=len-1;i>=0;i--)
  dseeds[i] = SUB2IND(xpts[i],ypts[i]);

setSeedPoints(dseeds, len);

::FREE(dseeds);
} 

#ifdef LDEBUG
#undef LDEBUG
#endif

#undef NROWS
#undef NCOLS

   /*------------------------ setSeedPoints ------------------------*/

/*-------------------------- setSeedRegions --------------------------*/
/*
   Defines the initial front to march out from.  Typically invoked with 
   a collection of points.  Can be invoked by specifying linear index
   locations or (x,y) image coordinate locations.

   void setSeedPoints( seeds, seed_length);
   void setSeedPoints( x coords, y coords, length);

   Inputs:
     seeds		- indices of seed points. 
     xpts, ypts		- coordinate version of seeds.
     slen		- number of the seed points

   Outputs: 		N/A.
   Requirements: 	N/A.
   Modifies:
     -the zero contour is defined according to the inputs.

   Memory Usage:
     -may allocate space for two integers if input is invalid.
     -allocates sufficient space to store the seed data.

     -for version taking xpts and ypts temporary space is allocated 
      to generate inputs that linear indices version of the setSeeds
      function accepts.
*/
void fastmarch::setSeedRegions(masktype* bwM, int mM, int mN)
{
register int i;

if (checkDimensions(mM, mN))
 {
  ERRORF("fastmarch::setImage -- dim Error",
         "Dimensions of seed matrix incompatible with expected dimensions.");
  return;
 }

for (i=area-1;i>=0;i--)
 {
  psi[i]  = 0.0;
  sign[i] = (bwM[i] >= 1)?ZERO:OUT;
 }

}  /*----------------------- setSeedRegions -----------------------*/

/*-------------------------- setKnownRegions -------------------------*/
/*
   Defines the initial front to march out from.  Typically invoked with 
   a collection of points, but in this case specified by a binary field
   over the domain.

   void setKnownRegions( masktype* bwM, int mM, int mN)

   Inputs:
     bwM		- the binary field (AKA, mask).
     mM, mN		- mask domain dimensions.

   Outputs: 		N/A.
   Requirements: 	N/A.
   Modifies:		N/A.
   Memory Usage:	N/A.
   Notes:
     as written, it may be a bit inefficient since two loops are used,
     one for distance initialization, and another for boundary
     initialization.  Both could be encapsulated in one loop.

	 in this case binary field means zero is empty and non-zero is a
	 known region.
*/
void fastmarch::setKnownRegions(masktype* bwM, int mM, int mN)
{
register int i;

if (checkDimensions(mM, mN))
 {
  ERRORF("fastmarch::setImage -- dim Error",
         "Dimensions of data incompatible with expected dimensions.");
  return;
 }

for (i=area-1;i>=0;i--)
 {
  psi[i]  = (bwM[i] == 0);
  sign[i] = (bwM[i] >= 1)?ZERO:OUT;
 }

compKnownBoundary();

}  /*----------------------- setKnownRegions ----------------------*/



/*------------------------------ setPhi ------------------------------*/
/*

*/

//TODO:  Add dimensions and check them?
void fastmarch::setPhi(number *nphi, int iM, int iN)
{

if (checkDimensions(iM, iN)) 
 {
  ERRORF("levelset::setImage -- dim Error",
         "Dimensions of data incompatible with expected dimensions.");
  return;
 }

if (ISNULL(phi))
 {
  phi = (number*)MALLOC(area*sizeof(number));

  if (ISNULL(phi))
   {
    ERRORF("fastmarch::init -- allocProb",
           "Can't allocate memory for fastmarch members.\n");
    return;
   }
 }

memcpy(phi, nphi, area*sizeof(number));

}


/*----------------------------- addtoPhi -----------------------------*/
void fastmarch::addtoPhi(number *aphi)
{
register int i;

for (i=area-1;i>=0;i--)
  phi[i] += aphi[i];

}

/*----------------------------- addtoPhi -----------------------------*/



/*----------------------------- scalePhi -----------------------------*/
void fastmarch::scalePhi(number *sphi)
{
register int i;

for (i=area-1;i>=0;i--)
  phi[i] *= sphi[i];

}

/*----------------------------- scalePhi -----------------------------*/



/*------------------------------ getPhi ------------------------------*/
/*

   Memory Usage:
     Allocates return space for Phi.  Must be freed manually
     outside of this function.
*/
number* fastmarch::getPhi(void)
{
number *retPhi;
register int i;

retPhi = (number *)MALLOC(area*sizeof(number));

if (!retPhi)
 {
  ERRORF("levelset::getPhi -- allocFail",
         "Could not allocate return space.");
 }
else
 memcpy(retPhi,phi,area*sizeof(number));

return retPhi;
}  /*--------------------------- getPhi --------------------------*/



#ifdef _MEX_
/*----------------------------- setParms -----------------------------*/
/*

*/

//#define LDEBUG 1	/* Local debug value. */

void fastmarch::setParms(mxArray *lsparms)
{
mxArray *fptr;		/* pointer to field element. */
double fval;		/* value of field element.   */
void *fdata;		/* pointer to field data, if array. */
int fM, fN;		/* field data dimensions, if array. */
number	*aphi;		/* additional potential field from workspace. */
int i;
int changePhi;


changePhi = 0;

matParms = lsparms;

//TODO: Implement a dirty field that determines if should run update or not.
//TODO:   Useful for callbacks, which may not need to update levelset
//TODO:   parameters, but just want to update evolution parameters.

/*--------- externally defined fields. ---------*/

fptr = mxGetField(lsparms,0,"phi");
if (fptr)
 {
  fdata = mxGetData(fptr);	/* additional potential component. */
  fM = (int)mxGetM(fptr);		/* num rows. */
  fN = (int)mxGetN(fptr);		/* num cols. */

  #if (defined(_DEBUG_) && defined(LDEBUG) && (LDEBUG))
  DEBUGF("  Overriding default potential field (AKA stopping function).\n");
  #endif
  setPhi((number *)fdata, fM, fN);
 }

fptr = mxGetField(lsparms,0,"sphi");
if (fptr)
 {
  fdata = mxGetData(fptr);		/* additional potential component. */
  fM = (int)mxGetM(fptr);		/* num rows. */
  fN = (int)mxGetN(fptr);		/* num cols. */

  #if (defined(_DEBUG_) && defined(LDEBUG) && (LDEBUG))
  DEBUGF("  Scaling internal potential using externally defined scaling.\n");
  #endif
  if ((M == fM) && (N == fN))
   {
    aphi = (number *)MALLOC(area*sizeof(number));
    for(i=area-1;i>=0;i--)
      aphi[i] = (number)(((double*)fdata)[i]);

    scalePhi(aphi);
    ::FREE(aphi);
   }
 }

fptr = mxGetField(lsparms,0,"aphi");
if (fptr)
 {
  fdata = mxGetData(fptr);		/* additional potential component. */
  fM = (int)mxGetM(fptr);		/* num rows. */
  fN = (int)mxGetN(fptr);		/* num cols. */

  #if (defined(_DEBUG_) && defined(LDEBUG) && (LDEBUG))
  DEBUGF("  Adding externally defined potential to internal potential.\n");
  #endif
  if ((M == fM) && (N == fN))
   {
    aphi = (number *)MALLOC(area*sizeof(number));
    for(i=area-1;i>=0;i--)
      aphi[i] = (number)(((double*)fdata)[i]);

    addtoPhi(aphi);
    ::FREE(aphi);
   }
 }

/*------------------ callback ------------------*/
#ifdef _MEX_
fptr = mxGetField(lsparms,0,"callback");
/*--Check for existence of the callback function?
    Doing so prevents seg faults when it doesn't exist.--*/
if (fptr)
 {
  fdata = (char*)mxArrayToString(fptr);		/* Name of callback function. */
  strncpy(callback, (char*)fdata, NAMELEN-1); 
  callback[NAMELEN-1] = 0x00;
  #if (defined(_DEBUG_) && defined(LDEBUG) && (LDEBUG))
  DEBUGF("  Using %s callback for evolution.\n", callback);
  //mexEvalString("exist('evolvePsep')");
  #endif
  ::FREE(fdata);
 }
#endif

}

#ifdef LDEBUG
#undef LDEBUG
#endif


    /*-------------------------- setParms --------------------------*/



/*---------------------------- mexPutData ----------------------------*/
/*
   
   Copies components of the levelset class to the Matlab workspace.

*/

void fastmarch::mexPutData(void)
{
mxArray *fptr;		/* pointer to field element. */
void *fdata;		/* pointer to data in fptr.  */

const char *field_names[] = {"psi","sign","phi","status","xdata","shocks"};
const int field_num = 6;
mwSize dims[2] = {1, 1};

mxArray *fmdata;

if (ISNULL(matVar))
 {
  matVar = mxCreateStructArray(1,dims,field_num,field_names);
  if (!matVar)
   {
    ERRORF("fastmarch::mexPutData -- createFail",
           "failed to create level set structure for workspace.");
    return;
   }
  if (isPersistent)
    ::ARRAYPERSIST(matVar);
 }

fmdata = matVar;

/*------------- copy psi to field. -------------*/
if (NOTNULL(psi))
 {
  fptr = mxGetField(fmdata,0,"psi");
  if (ISNULL(fptr))
   {
    fptr = mxCreateNumericMatrix(M, N, mxSINGLE_CLASS, mxREAL);
    if (ISNULL(fptr))
     {
      ERRORF("fastmarch::mexPutData -- allocFail",
             "failed to allocate variable for workspace.");
      return;
     }
   }
  fdata = mxGetPr(fptr);
  memcpy(fdata,psi,area*sizeof(number));
  mxSetField(fmdata, 0, "psi", fptr);
  if (isPersistent)
    ::ARRAYPERSIST(fptr);
 }

/*------------- copy sign to field. ------------*/
if (NOTNULL(sign))
 {
  fptr = mxGetField(fmdata,0,"sign");
  if (ISNULL(fptr))
   {
    fptr = mxCreateNumericMatrix(M, N, mxINT32_CLASS, mxREAL);
    if (ISNULL(fptr))
     {
      ERRORF("fastmarch::mexPutData -- allocFail",
             "failed to allocate variable for workspace.");
      return;
     }
   }
  fdata = mxGetPr(fptr);
  memcpy(fdata,sign,area*sizeof(signtype));
  mxSetField(fmdata, 0, "sign", fptr);
  if (isPersistent)
    ::ARRAYPERSIST(fptr);
 }

/*------------- copy phi to field. -------------*/
if (NOTNULL(phi))
 {
  fptr = mxGetField(fmdata,0,"phi");
  if (ISNULL(fptr))
   {
    fptr = mxCreateNumericMatrix(M, N, mxSINGLE_CLASS, mxREAL);
    if (ISNULL(fptr))
     {
      ERRORF("fastmarch::mexPutData -- allocFail",
             "failed to allocate variable for workspace.");
      return;
     }
   }
  fdata = mxGetPr(fptr);
  memcpy(fdata,phi,area*sizeof(number));
  mxSetField(fmdata, 0, "phi", fptr);
  if (isPersistent)
    ::ARRAYPERSIST(fptr);
 }

/*------------ copy status to field. -----------*/
if (NOTNULL(status))
 {
  fptr = mxGetField(fmdata,0,"status");
  if (ISNULL(fptr))
   {
    fptr = mxCreateNumericMatrix(M, N, mxINT32_CLASS, mxREAL);
    if (ISNULL(fptr))
     {
      ERRORF("fastmarch::mexPutData -- allocFail",
             "failed to allocate variable for workspace.");
      return;
     }
   }
  fdata = mxGetPr(fptr);
  memcpy(fdata,status,area*sizeof(fmstat));
  mxSetField(fmdata, 0, "status", fptr);
  if (isPersistent)
    ::ARRAYPERSIST(fptr);
 }

/*------------- copy xdat to field. ------------*/
if (NOTNULL(xdata))
 {
  fptr = mxGetField(fmdata,0,"xdata");
  if (ISNULL(fptr))
   {
    fptr = mxCreateNumericMatrix(M, N, mxINT32_CLASS, mxREAL);
    if (ISNULL(fptr))
     {
      ERRORF("fastmarch::mexPutData -- allocFail",
             "failed to allocate variable for workspace.");
      return;
     }
   }
  fdata = mxGetPr(fptr);
  memcpy(fdata,xdata,area*sizeof(fmstat));
  mxSetField(fmdata, 0, "xdata", fptr);
  if (isPersistent)
    ::ARRAYPERSIST(fptr);
 }

/*------------ copy shocks to field. -----------*/
if (NOTNULL(status))
 {
  fptr = mxGetField(fmdata,0,"shocks");
  if (ISNULL(fptr))
   {
    fptr = mxCreateNumericMatrix(M, N, mxINT32_CLASS, mxREAL);
    if (ISNULL(fptr))
     {
      ERRORF("fastmarch::mexPutData -- allocFail",
             "failed to allocate variable for workspace.");
      return;
     }
   }
  fdata = mxGetPr(fptr);
  memcpy(fdata,status,area*sizeof(fmstat));
  mxSetField(fmdata, 0, "shocks", fptr);
  if (isPersistent)
    ::ARRAYPERSIST(fptr);
 }



mexPutVariable("global",matName,fmdata);

}  /*------------------------- mexPutData -------------------------*/


/*---------------------------- mexPersist ----------------------------*/
/*
   Make memory persistent, if using mex implementation.  This is
   kind of annoying, but using Matlab memory management helps us
   out with memory problems.
*/
void fastmarch::mexPersist(void)
{

if (NOTNULL(psi))
  ::PERSIST(psi);
if (NOTNULL(sign))
  ::PERSIST(sign);
if (NOTNULL(phi))
  ::PERSIST(phi);
if (NOTNULL(status))
  ::PERSIST(status);
if (NOTNULL(shock))
  ::PERSIST(shock);
if (NOTNULL(xdata))
  ::PERSIST(xdata);
if (NOTNULL(seed))
  ::PERSIST(seed);

minheap.mexPersist();

isPersistent = true;

}  /*------------------------- mexPersist -------------------------*/
#endif


/*-------------------------- compZeroContour -------------------------*/
/*
   Find the zero level set contour.

   Returns:
     the length of the zero contour.

*/ 

#define POS OUT
#define NEG IN
#define NROWS M
#define NCOLS N

int fastmarch::compZeroContour(void)
{
short *marked;
int i, j, x, y;

marked = (short *)CALLOC(area,sizeof(short));
if ( !marked )
 {
  ERRORF("levelset::compZeroContour -- allocFail",
         "failed to allocate temporary storage.");
  return 0;
 }

seedlen = 0;

for(i=area-1;i>=0;i--)
 {
 
  if (sign[i] == ZERO)
   {
    seedlen++;
    marked[i] = 1;
   }
  else
   {
    x = XSUB(i);
    y = YSUB(i);

    if (x>0)
      if (sign[i] == POS && sign[i-DX] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }

    if (x<NCOLS-1)
      if (sign[i] == POS && sign[i+DX] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }

    if (y>0)
      if (sign[i] == POS && sign[i-DY] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }

    if (y<NROWS-1)
      if (sign[i] == POS && sign[i+DY] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }
   }
 }

if NOTNULL(seed)
  seed = (int *)REALLOC(seed,seedlen*sizeof(int));
else
  seed = (int *)MALLOC(seedlen*sizeof(int));

if (ISNULL(seed) && seedlen)
 {
  ERRORF("levelset::compZeroContour -- allocFail",
         "Unable to (re)allocate memory.");
  return -1;
 }

j = 0;
for(i=area-1;i>=0;i--)
  if (marked[i])
    seed[j++] = i;

::FREE(marked);

return seedlen;

}

int fastmarch::compZeroContour(number close)
{
register int i, j;
short *marked;

marked = (short *)CALLOC(area,sizeof(short));
if ( !marked )
 {
  ERRORF("levelset::compZeroContour -- allocFail",
         "failed to allocate temporary storage.");
  return 0;
 }

seedlen = 0;			/* Find close points. */
for (i=area-1;i>=0;i--)
  if (psi[i] <= close)
   {
    marked[i] = 1;
    seedlen++;
   }

if (NOTNULL(seed))
  seed = (int *)REALLOC(seed,seedlen*sizeof(int));
else
  seed = (int *)MALLOC(seedlen*sizeof(int));

if (ISNULL(seed) && seedlen)
 {
  ERRORF("levelset::compZeroContour -- allocFail",
         "Unable to (re)allocate memory.");
  return -1;
 }

j = 0;				/* Set close points to be the front. */
for(i=area-1;i>=0;i--)
  if (marked[i])
    seed[j++] = i;

::FREE(marked);

return seedlen;

#if defined(DEBUG) && defined(LDEBUG) && (LDEBUG)
if (LDEBUG)
  DEBUGF("The seed contour region within %5.2f is of area%d.\n", dist, zarea);
#endif

}

#undef POS
#undef NEG 
#undef NROWS
#undef NCOLS

   /*----------------------- compZeroContour ----------------------*/


/*-------------------------- compZeroRegion --------------------------*/
/*
   There are two methods for finding the zero region.

   (1) Find layer in and out of zero contour.  This is the ambiguous region
     within which the zero contour lies.  If a pixel is actually on the 
     zero contour, then just use the pixel and not the layer in and out 
     method.

   (2) Pass along an additional distance argument.  All values with
     this distance or smaller form the zero region.

   Memory Usage:
     Assumes that sufficient memory has been allocated to store the
     zero contour region in the argument zregion.  Typically this
     will be an chunk of memory with number elements equal to the 
     area of the level set.

*/

#define NROWS M
#define NCOLS N

#define LDEBUG 0

void fastmarch::compZeroRegion(int* zregion, int& zarea)
{
int i, x, y, marked;

zarea = 0;

for (i=area-1;i>=0;i--)
 {
  if (sign[i] == ZERO)
    zregion[zarea++] = i;
  else
   {
    marked = 0;
    x = XSUB(i);
    y = YSUB(i);

    if (x>0)
      if (sign[i]*sign[i-DX] < 0)
       {
        zregion[zarea++] = i;
	marked = 1;
       }

    if ((x<NCOLS-1) && !marked)
      if (sign[i]*sign[i+DX] < 0)
       {
        zregion[zarea++] = i;
	marked = 1;
       }

    if ((y>0) && !marked)
      if (sign[i]*sign[i-DY] < 0)
       {
        zregion[zarea++] = i;
	marked = 1;
       }

    if ((y<NROWS-1) && !marked)
      if (sign[i]*sign[i+DY] < 0)
       {
        zregion[zarea++] = i;
	marked = 1;
       }
   }      
 }

#ifdef DEBUG
if (LDEBUG)
  DEBUGF("Padded seed contour length is %d.\n",zarea);
#endif

}

void fastmarch::compZeroRegion(int* zregion, int& zarea, number dist)
{
register int i;

zarea = 0;

for (i=area-1;i>=0;i--)
  if (psi[i] <= dist)
    zregion[zarea++] = i;

#ifdef DEBUG
if (LDEBUG)
  DEBUGF("The seed contour region within %5.2f is of area%d.\n", dist, zarea);
#endif

}

#undef LDEBUG

   /*----------------------- compZeroRegion -----------------------*/


/*-------------------------- compKnownBoundary -------------------------*/
/*
   Find the boundary determined by the signed distance function.  Sets
   all inside points to KNOWN and all outside points to UNKNOWN during
   the procedure.  Also identifies the front region from which to propogate
   from.

   Returns:
     the length of the zero contour.

*/ 

#define POS OUT
#define NEG IN
#define NROWS M
#define NCOLS N

int fastmarch::compKnownBoundary(void)
{
short *marked;
int i, j, x, y;

marked = (short *)CALLOC(area,sizeof(short));
if ( !marked )
 {
  ERRORF("levelset::compZeroContour -- allocFail",
         "failed to allocate temporary storage.");
  return 0;
 }

seedlen = 0;

for(i=area-1;i>=0;i--)
 {
 
  if (sign[i] == ZERO)
   {
    seedlen++;
    marked[i] = 1;
   }
  else
   {
    x = XSUB(i);
    y = YSUB(i);

    if (x>0)
      if (sign[i] == POS && sign[i-DX] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }

    if (x<NCOLS-1)
      if (sign[i] == POS && sign[i+DX] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }

    if (y>0)
      if (sign[i] == POS && sign[i-DY] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }

    if (y<NROWS-1)
      if (sign[i] == POS && sign[i+DY] == NEG && !marked[i])
       {
        seedlen++;
        marked[i] = 1;
       }
   }
 }

if NOTNULL(seed)
  seed = (int *)REALLOC(seed,seedlen*sizeof(int));
else
  seed = (int *)MALLOC(seedlen*sizeof(int));

if (ISNULL(seed) && seedlen)
 {
  ERRORF("levelset::compZeroContour -- allocFail",
         "Unable to (re)allocate memory.");
  return -1;
 }

j = 0;
for(i=area-1;i>=0;i--)
  if (marked[i])
    seed[j++] = i;

::FREE(marked);

return seedlen;

}

#undef POS
#undef NEG 
#undef NROWS
#undef NCOLS

   /*----------------------- compKnownBoundary --------------------*/


/*======================== Distance Computation ========================*/
/* { */

/*---------------------------- initDist ----------------------------*/
/*
  Fill in the first layer of distances.  The remaining layers should 
  be filled in by fast marching or some other technique for propagating
  the distances out to the domain boundary.  Only 4-neighbors are 
  considered.  
  
  Prior to filling in the estimated distances, the entire grid status
  is reset to be FARAWAY (and at INFINITY if it is in the appropriate 
  direction). 

  Inputs:  		N/A.
  Outputs: 		N/A.
  Requirements: 	N/A.
  Modifies:
    -the distance map, psi, is reset and only the zero contour
     (and possibly its neighbors) reflect(s) the correct distances.

  Memory Usage:
    In the _SIMPLE_RI_ case, space is allocated for a copy of the
    distance function.  It is deallocated once no longer needed.
*/

#define NROWS M
#define NCOLS N

void fastmarch::initDist(void)
{
register int i, index;
int x, y, next;
heapType el;

/*------------- Initialize data. -------------*/

#if defined(_SIMPLE_RI_)
number *oldpsi;

oldpsi = (number*)MALLOC(seedlen*sizeof(number));/* allocate space (. */
if ( (seedlen > 0)  && ISNULL(oldpsi) )
 {
  ERRORF("levelset::initDist -- allocFail",
         "Could not allocate temporary space.");
  return;
 }
#endif

minheap.empty();				/* Prepare the heap. */
el.set(NUMBER(1), 0);

if (seedlen)					/* Done only if zero contour exists. */
 {
  #ifdef _SIMPLE_RI_
  for(i=seedlen-1; i>=0; i--)	/* Keep copy of dist. using pts. in zero. */
    oldpsi[i] = psi[seed[i]];
  #endif

  for(i=area-1; i>=0; i--)		/* Clean slate, all pts FAR at INFINITY, */
   {							/*   and no shocks.                      */
    status[i] = FAR;
    psi[i] = INFDIST;

    shock[i] = 0;
   }

  for(i=seedlen-1; i>=0; i--)	/* Start to fill in slate w/known data.  */
   {
    index = seed[i];
    status[index] = KNOWN;

    #if defined(_SIMPLE_RI_)
    psi[index] = oldpsi[i];		/* Keep old distance value. */
    #else
    psi[index] = 0;				/* Reset distance to be zero. */
    #endif
   }

 }

/*------- Any preparatory code goes here. ------*/

#if defined(_SIMPLE_RI_)

if (seedlen)
  ::FREE(oldpsi);				/* deallocate space ). */

#else

#endif


/*----------------- Main loop. -----------------*/

for(i=seedlen-1; i>=0; i--)
 {
  index = seed[i];

  #if defined(_SIMPLE_RI_)

  el.set(psi[index],index);	/* Setup zero region point heap element. */
  minheap.add(el);		/* Add element to heap. */

  #else

  x = XSUB(index);
  y = YSUB(index);

  next = index+DX;
  if ( (x<(NCOLS-1)) && (nogo[next] == 0) )
   {
    /* If FARAWAY, then this is first time we've seen this point.  Initialize
       as if it had only one neighbor.  If already marked as TRIAL, this point 
       touches at least two zero-level points, so let distance be sqrt(2)/2 
       instead of one.  Since sqrt(2)/2 is less than unity, the heap must 
       be adjusted. Of course, there are excceptions to this rule, but
       overall it works out OK (exceptions may be rare).
    */
    if (status[next] == FAR)
     {
      status[next] = TRIAL;
      psi[next] = NUMBER(1);
      el.setIndex(next);
      minheap.add(el);
     }
    else if (status[next] == TRIAL)
     {
      psi[next] = ROOT2OVER2;
      minheap.update(minheap.getId(next), (number)psi[next]);
     }
   }

  next = index-DX;
  if ( (x>0) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)
     {
      status[next] = TRIAL;
      psi[next] = NUMBER(1);
      el.setIndex(next);
      minheap.add(el);
     }
    else if (status[next] == TRIAL)
     {
      psi[next] = ROOT2OVER2;
      minheap.update(minheap.getId(next), (number)psi[next]);
     }
   } 


  next = index+DY;
  if ( (y<(NROWS-1)) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)
     {
      status[next] = TRIAL;
      psi[next] = NUMBER(1);
      el.setIndex(next);
      minheap.add(el);
     }
    else if (status[next] == TRIAL) 
     {
      psi[next] = ROOT2OVER2;
      minheap.update(minheap.getId(next), (number)psi[next]);
     }
   }

  next = index-DY;
  if ( (y>0) && (nogo[next] == 0) )
   {
    if (status[next] == FAR) 
     {
      status[next] = TRIAL;
      psi[next] = NUMBER(1);
      el.setIndex(next);
      minheap.add(el);
     }
    else if (status[next] == TRIAL)
     {
      psi[next] = ROOT2OVER2;
      minheap.update(minheap.getId(next), (number)psi[next]);
     }
   }

  #endif

 }

}

#undef NROWS
#undef NCOLS

   /*------------------------- initDist -------------------------*/



/*---------------------------- interpDist ----------------------------*/
/*
   Interpolates the location of the front using a piece-wise linear
   approximation.

   There are fives ways to have a 4-neighbor front crossing 
   (modulo rotation),

     (1) only one neighbor.

	  -------
          | |*| |   
	  -------
          | |+| |	- Just use interpolated distance.
	  -------
          | | | |
	  -------
	  

     (2) two contiguous neighbors

	  -------
          | |*| |   
	  -------
          |*|+| |	- d = \sqrt(  (1/s^2 + 1/t^2)^(-1) )
	  -------
          | | | |
	  -------
	  
     (3) three contiguous neighbors

	  -------
          | |*| |   
	  -------
          |*|+| |	- d = t*min(s_1,s_2) / sqrt(t^2 +  min^2(s_1,s_2) )
	  -------
          | |*| |
	  -------

     (4) two opposite neighbors

	  -------
          | |*| |   
	  -------
          | |+| |	- d = min(s_1, s_2)
	  -------
          | |*| |
	  -------

     (5) four neighbors

	  -------
          | |*| |   
	  -------                     min(s_1,s_2) min(t_1, t_2)
          |*|+|*|	- d =  --------------------------------------
	  -------              sqrt( min^2(s_1,s_2) + min^2(t_1,t_2))
          | |*| |
	  -------

   The rotation is determined by keeping track of who the 4-neighbors
   are using a boolean value system:

	  -------
          | |1| |   
	  -------
          |2|+|8|
	  -------
          | |4| |
	  -------
  

   Reference:
     Sethian, J.A., "Level Set Methods and Fast Marching Methods."

   Inputs:	N/A.
   Outputs:	N/A.
   Requires:	A valid distance function.

   Memory:	

   Notes:
     -The interpolated distances should be propogated over the
      entire domain by fast marching out the distances.
*/

#define NROWS M
#define NCOLS N

void fastmarch::interpDist(void)
{
register int i;
int x, y, xmax, ymax;
int   ncnt; 		/* Number of neighbors and       */
uint8 nbool;		/* their boolean interpretation. */

number dist_up, dist_dn,/* If near boundary, interpolated distances */ 
       dist_lf, dist_rt;/*   in the correct direction to boundary.  */
number psi_curr;	/* Current indexed value of level set.      */
signtype sgn_curr;	/* Current indexed sign of level set.       */
number *newpsi;		/* Pointer to new (allocated) distances.    */
number *psiptr, *newptr;/* Running pointers to distance functions.  */ 
signtype *sgnptr;	/* Running pointer to sign function.        */
fmstat *statptr;	/* Running pointer to fast march status pointer.  */
int8 *nogoptr;

number min_s, min_t;	/* For storing min. value comparisons of neighbors. */

heapType el;		/* Scratch heap element for dumping data to heap. */

newpsi = (number*)MALLOC(area*sizeof(number));       /* allocate space (. */
if (ISNULL(newpsi))
 {
  ERRORF("levelset::interpDist -- allocFail",
         "Could not allocate space for distance function.");
  return;
 }

xmax = NCOLS - 1;	/* Max index having neighbor to right. */
ymax = NROWS - 1;	/* Max index having neighbor below it. */

minheap.empty();	/* Prepare the heap. */

i = area - 1;		/* Set index to last element of array.        */
statptr = status + i;	/* Set all running pointers to last element   */
sgnptr  = sign + i;	/*   of their respective arrays.              */
newptr  = newpsi + i;
psiptr  = psi + i;
nogoptr = nogo + i;


for (; i>=0; i--)
 {
  sgn_curr = *sgnptr;			/* Current value of the sign. */
  if (sgn_curr == ZERO)
   {
    (*newptr)  = (number)0;
    (*statptr) = KNOWN;

    el.set( (number)0, i);		/* Update heap element data. */
    minheap.add(el);			/* Copy element to the heap. */
   }
  else
   {
    ncnt = 0;				/* Reset neighbor count. */
    nbool = 0x00;			/* Reset neighbor boolean. */
    psi_curr = *psiptr;			/* Current distance value. */

    x = XSUB(i);			/* Get coordinate location of index.  */
    y = YSUB(i);

    if ( (y > 0) && (*(nogoptr-DY) == 0) && (sgn_curr*(*(sgnptr-DY)) < 0) )
     {
      ncnt++;
      nbool += 0x01;			/* Add an up neighbor. */
      dist_up = psi_curr/(psi_curr + *(psiptr-DY));
     }

    if ( (x > 0) && (*(nogoptr-DX) == 0) && (sgn_curr*(*(sgnptr-DX)) < 0) )
     {
      ncnt++;
      nbool += 0x02;			/* Add a left neighbor. */
      dist_lf = psi_curr/(psi_curr + *(psiptr-DX));
     }

    if ( (y < ymax) && (*(nogoptr+DY) == 0) && (sgn_curr*(*(sgnptr+DY)) < 0) )
     {
      ncnt++;
      nbool += 0x04;			/* Add a down neighbor. */
      dist_dn = psi_curr/(psi_curr + *(psiptr+DY));
     }

    if ( (x < xmax) && (*(nogoptr+DX) == 0) && (sgn_curr*(*(sgnptr+DX)) < 0) )
     {
      ncnt++;
      nbool += 0x08;			/* Add a right neighbor. */
      dist_rt = psi_curr/(psi_curr + *(psiptr+DX));
     }


    if (ncnt > 0)
     {
      if (ncnt == 1)			/* Only one neighbor. */
       {
        switch(nbool)
         {
          case 0x01:		/* up. */ 
		    (*newptr) = dist_up; 
		    break;
		  case 0x02:		/* left. */ 
		    (*newptr) = dist_lf; 
			break;
	      case 0x04:		/* down. */ 
		    (*newptr) = dist_dn; 
			break;
	      case 0x08:		/* right. */ 
		    (*newptr) = dist_rt; 
			break;
         }
       }
      else if (ncnt == 2)		/* Have two neighbors, two cases. */
       {
        switch (nbool) 
         {
          case 0x03:		/* up, left. */
	        (*newptr) = (dist_up * dist_lf ) /
	                      SQRT( (dist_up*dist_up) + (dist_lf*dist_lf) );
	        break; 
		  case 0x06:		/* down, left. */ 
		    (*newptr) = (dist_dn * dist_lf) /
	                      SQRT( (dist_dn*dist_dn) + (dist_lf*dist_lf) );
	        break;
	      case 0x09:		/* up, right. */
	        (*newptr) = (dist_up * dist_rt) /
	                      SQRT( (dist_up*dist_up) + (dist_rt*dist_rt) );
	        break;
	      case 0x0C:		/* down, right. */
	        (*newptr) = (dist_dn * dist_rt) / 
	                      SQRT( (dist_dn*dist_dn) + (dist_rt*dist_rt) );
	        break;
	      case 0x05:		/* up, down. */
	        (*newptr) = MIN(dist_up, dist_dn);
	        break;
	      case 0x0A:		/* left, right. */
	        (*newptr) = MIN(dist_lf, dist_rt);
	        break;
         }
       }
      else if (ncnt == 3)
       {
        switch(nbool)
         {
          case 0x07:		/* up, left, down. */
            min_s = MIN(dist_up, dist_dn);
            (*newptr) = (min_s * dist_lf) / 
                                    SQRT( (dist_lf*dist_lf) + (min_s*min_s) );
            break;
          case 0x0B:		/* right, up, left. */
            min_s = MIN(dist_lf, dist_rt);
            (*newptr) = (min_s * dist_up) /
                                    SQRT( (dist_up*dist_up) + (min_s*min_s) );
            break;
          case 0x0D:		/* down, right, up. */
            min_s = MIN(dist_up, dist_dn);
            (*newptr) = (min_s * dist_rt) /
                                    SQRT( (dist_rt*dist_rt) + (min_s*min_s) );
            break;
          case 0x0E:		/* left, down, right. */
            min_s = MIN(dist_lf, dist_rt);
            (*newptr) = (min_s * dist_dn) / 
                                    SQRT( (dist_dn*dist_dn) + (min_s*min_s) );
            break;
         }
       }
      else if (ncnt == 4)	/* all are neighbors. */
       {
        min_s = MIN(dist_lf, dist_rt);
		min_t = MIN(dist_up, dist_dn); 
		(*newptr) = (min_s * min_t) / SQRT( (min_s*min_s) + (min_t*min_t) );
       }

      (*statptr) = KNOWN;		/* Interpolated distances are known. */
      el.set( (*newptr), i);		/* Update heap element data. */
      minheap.add(el);			/* Copy element to the heap. */
      #if defined(LDEBUG) && (LDEBUG >= 2)
      DEBUGF("Added distance %5.3f with sign of %d.\n", (*newptr), sign[i]);
      #endif
     }
    else
     {
      (*statptr) = FAR;			/* If not near front, then far away. */
      (*newptr)  = INFDIST;		/* Set distance to infinity. */
     }
   }

  statptr--;				/* Update running pointers. */
  sgnptr--;
  newptr--;
  psiptr--;
  nogoptr--;
 }

::FREE(psi);       				/* De-allocate old space ). */
psi = newpsi;					/* Keep the new distances.  */

}

#undef NROWS
#undef NCOLS

#ifdef LDEBUG
#undef LDEBUG
#endif

   /*------------------------- interpDist -------------------------*/



/*----------------------------- solveDist ----------------------------*/
/* 
   Solve for T at index using given speed function. Returns the new
   time value for the fast marching Boundary Value problem.  This
   function is specific to the fast march problem which aims to solve
   for the distances.  Only 4-neighbors are considered.

   The algorithm used to solve for the arrival time is first-order
   close. 
   
   Reference:
     Sethian, J.A., "Level Set Methods and Fast Marching Methods."

   Inputs:
     -index to the point to solve for.
     -speed of the propogating front.

   Outputs:
     -estimated arrival time of the point.

   Requirements: N/A.
   Memory Usage: N/A.
*/

#define NROWS M
#define NCOLS N

number fastmarch::solveDist(int index, number speed)
{
number T, Tup, Tdn, Tlf, Trt;		/* neighboring T values.  */
number Dup, Ddn, Dlf, Drt;		/* diffs in 4 directions. */
number A, B, C;				/* to solve quad. eqn.    */
float Disc;				/* discriminant 	  */
int x, y;


T = psi[index];			/* Grab current arrival time. */

x = XSUB(index);		/* Location of point for boundary check. */
y = YSUB(index);

if (x>0)			/* Compute neighboring arrival times. */
  Tlf = psi[index-DX];
else
  Tlf = psi[index];

if (x<(NCOLS-1))
  Trt = psi[index+DX];
else
  Trt = psi[index];

if (y>0)
  Tup = psi[index-DY];
else
  Tup = psi[index];

if (y<(NROWS-1))
  Tdn = psi[index+DY];
else
  Tdn = psi[index];

Dup = T-Tup; 			/* Compute discrete derivatives. */
Ddn = Tdn-T;
Dlf = T-Tlf;
Drt = Trt-T;

A = NUMBER(0); 			/* Setup quadratic equation. */
B = NUMBER(0);
C = -INV(TIMES(speed,speed));

if ((Dlf>0) || (Drt<0))		/* Check diffs. for proper upwind solution. */
 {							/*   also finds winner in case of shock.    */
  A++;

  if (Dlf > -Drt)
   {
    B -= TIM2(Tlf);
    C += TIMES(Tlf,Tlf);
   }
  else
   {
    B -= TIM2(Trt);
    C += TIMES(Trt,Trt);
   }
 }

if ((Dup>0) || (Ddn < 0))
 {
  A++;

  if (Dup > -Ddn)
   {
    B -= TIM2(Tup);
    C += TIMES(Tup,Tup);
   }
  else
   {
    B -= TIM2(Tdn);
    C += TIMES(Tdn,Tdn);
   }
 }

/* If all four differences were zero, A will be zero, and we cannot
   do the computation, but neither can we return the original T 
   value.  Instead we will take the shortest arrival time and compute
   an approximate arrival time based on that data point alone.
   In the quadratic solution for T, the plus sign is correct due to 
   upwind nature of problem.   Should rarely or never happen (this may
   be a legacy conditional that may be removable).
*/
if (A!=0)
 {
  Disc = TIMES(B,B) - POW2(TIMES(A,C),2);
  if (Disc > 0)
   {
    T = DIVIDE( -B+NUMBER(sqrt(Disc)) , TIM2(A) );
   }
  else	
   {
    if (T > Tup)
      T = Tup;
    if (T > Tdn)
      T = Tdn;
    if (T > Trt)
      T = Trt;
    if (T > Tlf)
      T = Tlf;
    T = T+INV(speed);
   }
 }
else
 {
  if (T > Tup)
    T = Tup;
  if (T > Tdn)
    T = Tdn;
  if (T > Trt)
    T = Trt;
  if (T > Tlf)
    T = Tlf;
  T = T+INV(speed);
 }

return T;
}

#undef NROWS
#undef NCOLS

/*------------------------- solveDist -------------------------*/

/*---------------------------- solveDistL1 ---------------------------*/
/* 
   Solve for T at index using given speed function. Returns the new
   time value for the fast marching Boundary Value problem.  This
   function is specific to the L1 fast march problem which aims to solve
   for the L1 distances.  Naturally, only 4-neighbors are considered.

   Since it uses the heap, it is a region growing method, but in reality
   it could be implemented as a Chamfer transform, which would make it
   operate in something like O(4*area) speed.

   Inputs:
     -index to the point to solve for.
     -speed of the propogating front.

   Outputs:
     -estimated arrival time of the point.

   Requirements: N/A.
   Memory Usage: N/A.
*/

#define NROWS M
#define NCOLS N

number fastmarch::solveDistL1(int index, number speed)
{
number T, Tup, Tdn, Tlf, Trt;		/* neighboring T values.   */
int x, y;


T = psi[index];			/* Grab current arrival time. */

x = XSUB(index);		/* Location of point for boundary check. */
y = YSUB(index);

if (x>0)			/* Compute neighboring arrival times. */
  Tlf = psi[index-DX];
else
  Tlf = psi[index];

if (x<(NCOLS-1))
  Trt = psi[index+DX];
else
  Trt = psi[index];

if (y>0)
  Tup = psi[index-DY];
else
  Tup = psi[index];

if (y<(NROWS-1))
  Tdn = psi[index+DY];
else
  Tdn = psi[index];

if (T > Tup)
  T = Tup;
if (T > Tdn)
  T = Tdn;
if (T > Trt)
  T = Trt;
if (T > Tlf)
  T = Tlf;

T = T+INV(speed);

return T;
}

#undef NROWS
#undef NCOLS

   /*------------------------- solveDistL1 ------------------------*/

/*---------------------------- solveDistN8 ---------------------------*/
/* 
   Solve for T at index using given speed function. Returns the new
   time value for the fast marching Boundary Value problem.  The algorithm
   looks at the 8-neighbors and picks the neighbor with the minimum value.
   This Chamfer distance method is clearly an approximation only.

   Method 1 gives a maximum error of:	7.61%	(Octoganal Euclidean metric)
   Method 2 gives a maxumum error of:	3.96%	(Butt and Maragos, 1998)

   Can be made more efficient by using a double-scan approach.  There is no
   need for a heap.  It was just really fast to modify this way.  It
   should operate in something like O(2*m*area) where m is the number of
   neighboring cells looked at in the half-stencil.

   Inputs:
     -index to the point to solve for.
     -speed of the propogating front.

   Outputs:
     -estimated arrival time of the point.

   Requirements: N/A.
   Memory Usage: N/A.
*/

#define NROWS M
#define NCOLS N

#ifndef _CHAMFDIST_
#define _CHAMFDIST_ 2
#endif

#if (_CHAMFDIST_ == 1)
#define CHAMF_A		1
#define CHAMF_B		ROOT2
#elif (_CHAMFDIST_ == 2)
#define CHAMF_A		0.96194
#define CHAMF_B		1.36039
#endif


number fastmarch::solveDistN8(int index, number speed)
{
number Tup, Tul, Tur,
       T  , Tlf, Trt,			/* neighboring T values.   */
       Tdn, Tdl, Tdr;
int x, y;
number fact;


T = psi[index];			/* Grab current arrival time. */

x = XSUB(index);		/* Location of point for boundary check. */
y = YSUB(index);

if (x>0)			/* Compute neighboring arrival times. */
 {
  Tlf = psi[index-DX];

  if (y>0)
    Tul = psi[index-DX-DY];
  else
    Tul = psi[index];

  if (y<(NCOLS-1))
    Tdl = psi[index-DX+DY];
  else
    Tdl = psi[index];
 }
else
 {
  Tlf = psi[index];
  Tul = psi[index];
  Tdl = psi[index];
 }

if (x<(NCOLS-1))
 {
  Trt = psi[index+DX];

  if (y>0)
    Tur = psi[index+DX-DY];
  else
    Tur = psi[index];

  if (y<(NCOLS-1))
    Tdr = psi[index+DX+DY];
  else
    Tdr = psi[index];
 }
else
 {
  Trt = psi[index];
  Tur = psi[index];
  Tdr = psi[index];
 }

if (y>0)
  Tup = psi[index-DY];
else
  Tup = psi[index];

if (y<(NROWS-1))
  Tdn = psi[index+DY];
else
  Tdn = psi[index];

fact = CHAMF_A;
if (T > Tup)
  T = Tup;
if (T > Tdn)
  T = Tdn;
if (T > Trt)
  T = Trt;
if (T > Tlf)
  T = Tlf;

if (T > Tur)
 {
  T = Tur;
  fact = CHAMF_B;
 }
if (T > Tul)
 {
  T = Tul;
  fact = CHAMF_B;
 }
if (T > Tdl)
 {
  T = Tdl;
  fact = CHAMF_B;
 }
if (T > Tdr)
 {
  T = Tdr;
  fact = CHAMF_B;
 }

T = T+fact*INV(speed);

return T;
}

#undef NROWS
#undef NCOLS

   /*------------------------- solveDistN8 ------------------------*/


/*------------------------- setDistance -------------------------*/
/* sets the psi and the sign field from a signed distance function */
//TODO: what's the purpose.

void fastmarch::setDistance( number *sdist ) 
{
int i;
number psi_i;

for ( i=area-1; i>=0; i-- ) 
 {
  psi_i = sdist[i];

  if ( psi_i>0 ) 
   {
    psi[i] = psi_i;
    sign[i] = OUT;
   } 
  else if ( psi_i<0 ) 
   {
    psi[i] = -psi_i;
    sign[i] = IN;
   } 
  else 
   {
    psi[i] = 0.0;
    sign[i] = ZERO;
   }
  }

//TODO: Is this necessary.  May be if this function is even called.
//TODO: But what if we want to use alternative level set values.
//TODO: May then want some kind distance set with offset (for other levels)
//fastmarch::compZeroContour();

}  /*---------------------- setDistance ----------------------*/

/*----------------------------- distmarch -----------------------------*/
/*
  This function fast marches out an initialized distance map.  While the
  heap is not empty, there are points whose corresponding distances
  need to be updated.  Once the heap is empty, there is no more work
  to do and the distance map should be bonafide.

  Inputs: N/A.
  Outputs: N/A.
  Requirements:
    -the distance map should have been prepared initially be initDist.

  Memory Usage: N/A.

*/
#define NROWS M
#define NCOLS N

#define LDEBUG  0	/*  Local debug value.  Adjust to needs. */
#ifndef SOLVEDIST
#define SOLVEDIST	solveDist
#endif

void fastmarch::distmarch(void)	 
{
register number T;
register int next;
int mindex, x, y;
unsigned char testH, testV, testD1, testD2;
heapType el;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

//int iloop = minheap.numel();
//while (iloop-- > 0) 
while (minheap.pull(el))
 { 
  mindex = el.getIndex();
  x = XSUB(mindex);
  y = YSUB(mindex);

  testH = 0x00;
  testV = 0x00;
  testD1 = 0x00;
  testD2 = 0x00;


  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed (%d, %d) value %f.\n", ++cnt, x, y, 
           el.getTime());
  #endif

  next = mindex-DX;
  if (x>0)
   {
    if (nogo[next] == 0)
	 {
      if (status[next] == FAR)			/* add to heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
        status[next] = TRIAL;
        psi[next] = T;
        el.set(T,next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 2)
          DEBUGF("\t1--Added %d (%d, %d) with value %f.\n", next, x-1, y, T);
        #endif
       }
      else if (status[next] == TRIAL)		/* update heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t1--Updated %d (%d,%d) from %f to %f.\n",next, x-1, y, 
  	       psi[next], T);
        #endif
  
        psi[next] = T;
        minheap.update(minheap.getId(next),T);
       }
  
      if (status[mindex] != KNOWN)
       {
        if (psi[next] == psi[mindex])
          testH = 1;
        else if (psi[next] < psi[mindex])
          testH = 2;
       }
     }

    if (y>0)
     {
      next-=DY;
	  if (nogo[next] == 0)
	   {
        if (psi[next] == psi[mindex])
          testD1 = 1;
        else if (psi[next] < psi[mindex])
          testD1 = 2;
	   }
     }
    if (y<(NROWS-1))
     {
	  if (nogo[next] == 0)
	   {
        next+=(DY+DY);
        if (psi[next] == psi[mindex])
          testD2 = 1;
        else if (psi[next] < psi[mindex])
        testD2 = 2;
	   }
     }
   }


  next = mindex+DX;
  if (x<(NCOLS-1))
   {
    if (nogo[next] == 0)
	 {
      if (status[next] == FAR)			/* add to heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
        status[next] = TRIAL;
        psi[next] = T;
        el.set(T,next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 2)
          DEBUGF("\t2--Added %d (%d, %d) with value %f.\n",next, x+1, y, T);
        #endif
  
       }
      else if (status[next] == TRIAL)		/* update heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t2--Updated %d (%d,%d) from %f to %f.\n",next, x+1, y, 
  	       psi[next], T);
        #endif
  
        psi[next] = T;
        minheap.update(minheap.getId(next),T);
       }
  
      if (testH)
       {
        if ( (testH == 1) && (psi[next] < psi[mindex]) )
          testH = 5;
        else if ( (testH == 2) && (psi[next] <= psi[mindex]) )
          testH = 5;
       }
     }

    if (y<(NROWS-1))
     {
      next+=DY;
	  if (nogo[next] == 0)
	   {
        if ( (testD1 == 1) && (psi[next] < psi[mindex]) )
          testD1 = 5;
        else if ( (testD1 == 2) && (psi[next] <= psi[mindex]) )
          testD1 = 5;
	   }
     }
    if (y>0)
     {
      next-=(DY+DY);
	  if (nogo[next] == 0)
	   {
        if ( (testD2 == 1) && (psi[next] == psi[mindex]) )
          testD2 = 5;
        else if ( (testD2 == 2) && (psi[next] < psi[mindex]) )
          testD2 = 5;
	   }
     }
   }

  next = mindex-DY;
  if ( (y>0) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",next, x, y-1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Updated %d (%d,%d) from %f to %f.\n",next, x, y-1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }

    if (status[mindex] != KNOWN)
     {
      if (psi[next] == psi[mindex])
        testV = 1;
      else if (psi[next] < psi[mindex])
        testV = 2;
     }
   }

  next = mindex+DY;
  if ( (y<(NROWS-1)) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", next, x, y+1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Updated %d (%d,%d) from %f to %f.\n",next, x, y+1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }

    if (testV)
     {
      if ( (testV == 1) && (psi[next] < psi[mindex]) )
        testV = 5;
      else if ( (testV == 2) && (psi[next] <= psi[mindex]) )
        testV = 5;
     }
   }

  shock[mindex] = (testH == 5) + ((testV == 5)<<1) 
                    + ((testD1 ==5)<<2) + ((testD2==5)<<3);
  status[mindex] = KNOWN;
 } 

}

#undef LDEBUG

#undef NROWS
#undef NCOLS

   /*------------------------- distmarch -------------------------*/

/*---------------------------- fieldmarch ----------------------------*/
/*
  This function fast marches out an initialized distance map, where the
  traversal speed is given by a scalar field (speed varies spatially).  
  While the heap is not empty, there are points whose corresponding 
  distances need to be updated.  Once the heap is empty, there is no 
  more work to do and the traversal time map should be bonafide.

  Inputs: N/A.
  Outputs: N/A.
  Requirements:
    -the traversal speed map should have been prepared initially be initDist.

  Memory Usage: N/A.


  TODO: Does not generate a shock map like distmarch does.  Needs coding.

*/
#define NROWS M
#define NCOLS N

#define LDEBUG 0	/*  Local debug value.  Adjust to needs. */

void fastmarch::fieldmarch(number *speed)
{
register number T;
register int next;
int mindex, x, y;
heapType el;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

while (minheap.pull(el))
 { 
  mindex = el.getIndex();
  status[mindex] = KNOWN;
  x = XSUB(mindex);
  y = YSUB(mindex);

  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed (%d, %d) value %f.\n", ++cnt, x, y, 
           el.getTime());
  #endif

  next = mindex-DX;
  if ( (x>0) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = SOLVEDIST(next, speed[mindex]);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t1--Added %d (%d, %d) with value %f.\n", next, x-1, y, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Updated %d (%d,%d) from %f to %f.\n",next, x-1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }

  next = mindex+DX;
  if ( (x<(NCOLS-1)) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = SOLVEDIST(next, speed[mindex]);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t2--Added %d (%d, %d) with value %f.\n",next, x+1, y, T);
      #endif

     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Updated %d (%d,%d) from %f to %f.\n",next, x+1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }

  next = mindex-DY;
  if ( (y>0) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = SOLVEDIST(next, speed[mindex]);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",next, x, y-1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Updated %d (%d,%d) from %f to %f.\n",next, x, y-1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }

  next = mindex+DY;
  if ( (y<(NROWS-1)) && (nogo[next] == 0) )
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = SOLVEDIST(next, speed[mindex]);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", next, x, y+1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Updated %d (%d,%d) from %f to %f.\n",next, x, y+1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }
 } 

}

#undef LDEBUG

#undef NROWS
#undef NCOLS

   /*------------------------- fieldmarch -------------------------*/

/* } */
/*======================== Distance Computation ========================*/

/*========================== Label Propogation =========================*/
/* { */

/*-------------------------- setLabelPoints --------------------------*/
/*
   Defines the points in space that are associated with a label.  It
   is further presumed that these points will also form the initial front 
   points/regions to march out from.  Each point is assumed to be an
   individual label.  To give multiple points a common label, use
   setLabelRegions instead (which requires a label map).

   void setLabelPoints( seeds, seed_length);
   void setLabelPoints( x coords, y coords, length);

   Inputs:
     dlabels		- the label indices.
	 dlen			- the number of label indices.

   Outputs: 		N/A.

   Requirements: 	N/A.

   Modifies:
     -the seed points are defined according to the label points.
     -the distance map is wiped clean.

   Memory Usage:
     -allocates sufficient space to store the label data, if not done.

     -for version taking xpts and ypts temporary space is allocated 
      to generate inputs that linear indices version of the function 
	  accepts.
*/

#define NROWS M
#define NCOLS N

//#define LDEBUG 2

void fastmarch::setLabelPoints(int *dlabels, int dlen)
{
int *newptr;
int i, j;

newptr = NULL;


/*--(1) First deal with the contour seeds. --*/
if (seedlen > 0)					/* if zero contour exists, wipe it.     */
 {
  seedlen = 0;
  ::FREE(seed);
 }

if (dlen <= 0)						/* if passed length is nonsensical ...  */
 {
  seed = (int *)MALLOC(1*sizeof(int));
  if (ISNULL(seed))
   {
    ERRORF("fastmarch::setLabelPoints::allocFail",
           "Could not allocate seed memory.");
    return;
   }
  *seed = area>>1; 					/* set seed point to center of image.   */
  seedlen = 1;						
 }
else
 {
  seedlen = dlen;					/* otherwise, copy desired into actual. */

  seed = (int *) MALLOC(seedlen*sizeof(int));
  if (ISNULL(seed))
   {
    ERRORF("fastmarch::setLabelPoints::allocFail",
           "Unable to allocate memory for Zero Contour.");
    return;
   }

  memcpy(seed,dlabels,dlen*sizeof(int));
 }


/*--(2) Then deal with the label map allocation.  --*/
if ( ISNULL(labels) )
  labels = (labeltype *) CALLOC(area, sizeof(labeltype));

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::setLabelPoints -- allocProb",
         "Can't allocate memory for fastmarch member labels.\n");
  return;
 }


/*--(3) Lastly, wipe the slate clean and initialize. --*/
for (i=seedlen-1;i>=0;i--)			/* for each element of zero contour ... */
 {
  j = seed[i];						/* get index value,                 */
  psi[j]  = 0;						/* set distance to zero,            */
  sign[j]  = ZERO;
  status[j] = KNOWN;				/* fast march status known.         */
  labels[j] = i+1;					/* set label according to ordering. */
 }

}

void fastmarch::setLabelPoints(int *xpts,int *ypts, int len)
{
int *dlabels, i;

dlabels = (int *)MALLOC(sizeof(int)*len);

if (ISNULL(dlabels))
 {
  ERRORF("fastmarch::setLabelPoints::allocFail",
         "Unable to allocate memory for label index list.");
  return;
 }

for (i=len-1;i>=0;i--)
  dlabels[i] = SUB2IND(xpts[i],ypts[i]);

setLabelPoints(dlabels, len);

::FREE(dlabels);
}

#ifdef LDEBUG
#undef LDEBUG
#endif

#undef NROWS
#undef NCOLS

   /*----------------------- setLabelPoints -----------------------*/

/*-------------------------- setLabelRegions -------------------------*/
/*
   Defines the initial front regions to march out from and the labels
   associated to the regions.  Labels are presumed to be greater than
   zero.

   Inputs:
     lmap		- the label map.
     mM, mN		- the label map dimensions.

   Outputs: 		N/A.
   Requirements: 	N/A.
   Modifies:
     -the seed region/contour is defined according to the labelled 
      regions.
     -the distance map is wiped clean.

   Memory Usage:
     -allocates sufficient space to store the label data.

     -for version taking xpts and ypts temporary space is allocated 
      to generate inputs that linear indices version of the setSeeds
      function accepts.
*/
void fastmarch::setLabelRegions(labeltype* nlabel, int mM, int mN)
{
register int i;

if (checkDimensions(mM, mN))
 {
  ERRORF("fastmarch::setImage -- dim Error",
         "Dimensions of data incompatible with expected dimensions.");
  return;
 }

if ( ISNULL(labels) )
  labels = (labeltype *) MALLOC(area*sizeof(labeltype));

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::setLabelRegions -- allocProb",
         "Can't allocate memory for fastmarch member labels.\n");
  return;
 }

for (i=area-1;i>=0;i--)
 {
  psi[i]   = 0.0;
  sign[i]  = (nlabel[i] > 0)?ZERO:OUT;
  labels[i] = nlabel[i];
 }

}  /*----------------------- setLabelRegions ----------------------*/



/*---------------------------- getLabels -----------------------------*/
/*

   Memory Usage:
     Allocates return space for the label map.  Must be freed manually 
     outside of this function.
*/
labeltype* fastmarch::getLabels(void)
{
labeltype *retlabels;

retlabels = (labeltype *)MALLOC(area*sizeof(labeltype));

if (ISNULL(retlabels))
 {
  ERRORF("levelset::getLabels -- allocFail",
         "Could not allocate return space for labels.");
 }
else
  memcpy(retlabels, labels, area*sizeof(labeltype));

return retlabels;
}

void fastmarch::getLabels(labeltype *rlabels)
{
memcpy(rlabels, labels, area*sizeof(labeltype));
}

   /*------------------------- getLabels --------------------------*/


/*---------------------------- compLabels ----------------------------*/
/*
*/
void fastmarch::compLabels()
{
relabel();
}  /*------------------------- compLabels -------------------------*/

/*------------------------ compLabelsCostField -----------------------*/
/*
    Memory:
      Allocates space to convert desired cost field into an effective
      cost field for the fast marching function.  The presumed
      relationship is inverse.  If other monotonic relationships are
      desired for conversion between cost and speed, then conversion should
      be done beforehand with compdistSpeedField as the proper marching
      function.
*/

void fastmarch::compLabelsCostField( number *cost )
{
int i;
number *speed;

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::relabel -- missingInit",
         "Need to first initialize/define the initial label map.\n");
  return;
 }

speed = (number *)MALLOC(area*sizeof(number));
if (ISNULL(speed))
 {
  ERRORF("levelset::compdistCostField -- allocFail",
         "Could not allocate return space.");
 }
else
 {
  for (i=area-1;i>=0;i--)
    speed[i] = 1/cost[i];
 }

compZeroContour(); 
initDist();
labelfieldmarch(speed);
}  
   /*--------------------- compLabelsCostField --------------------*/


/*----------------------- compLabelsSpeedField -----------------------*/
/*
*/

void fastmarch::compLabelsSpeedField( number *speed )
{
int i;

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::relabel -- missingInit",
         "Need to first initialize/define the initial label map.\n");
  return;
 }

compZeroContour(); 
initDist();
labelfieldmarch(speed);
}  
   /*-------------------- compLabelsSpeedField --------------------*/

/*------------------------- compLabelsByCost -------------------------*/
/*
    Memory:
      Allocates space to convert desired cost field into an effective
      cost field for the fast marching function.  The presumed
      relationship is inverse.  If other monotonic relationships are
      desired for conversion between cost and speed, then conversion should
      be done beforehand with compdistSpeedField as the proper marching
      function.
*/

void fastmarch::compLabelsByCost( number *cost , int numlabels)
{
int i;
number *speed;

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::relabel -- missingInit",
         "Need to first initialize/define the initial label map.\n");
  return;
 }

speed = (number *)MALLOC(area*sizeof(number));
if (ISNULL(speed))
 {
  ERRORF("levelset::compdistCostField -- allocFail",
         "Could not allocate return space.");
 }
else
 {
  for (i=numlabels-1;i>=0;i--)
   {
    speed[i] = 1/cost[i];
	printf("Cost = %f, Speed = %f.\n", cost[i], speed[i]);
   }
 }

compZeroContour(); 
initDist();
labelspeedmarch(speed);
}  
   /*---------------------- compLabelsByCost ----------------------*/


/*------------------------- compLabelsBySpeed ------------------------*/
/*
*/

void fastmarch::compLabelsBySpeed( number *speed , int numlabels)
{
int i;

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::relabel -- missingInit",
         "Need to first initialize/define the initial label map.\n");
  return;
 }

compZeroContour(); 
initDist();
labelspeedmarch(speed);
}  
   /*---------------------- compLabelsBySpeed ---------------------*/


/*---------------------------- labelmarch ----------------------------*/
/*
   Fast march out the label terms to identify the regions associated to 
   the different "sources."

   Inputs: N/A.

   Outputs: N/A.

   Requirements:
     -the heap and initial labels must be initialized.

   Modifies:
     -the extension data is propogated over the entire domain.
     -the levelset heap should be empty.

   Memory Usage: N/A.

*/
#define NROWS M
#define NCOLS N

#define LDEBUG  0	/*  Local debug value.  Adjust to needs. */

void fastmarch::labelmarch(void)
{
register number T;
register int next;
int mindex, x, y;
labeltype lab_lf, lab_rt, lab_up, lab_dn;
labeltype lab_ul, lab_dl, lab_ur, lab_dr;
heapType el;

number bpsi;
labeltype blabel;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

while (minheap.pull(el))
 { 
  mindex = el.getIndex();
  x = XSUB(mindex);
  y = YSUB(mindex);

  blabel = 0;
  lab_lf = 0;
  lab_rt = 0;
  lab_up = 0;
  lab_dn = 0;
  lab_ul = 0;
  lab_dl = 0;
  lab_ur = 0;
  lab_dr = 0;

  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed (%d, %d) value %f.\n", ++cnt, x, y, 
           el.getTime());
  #endif

  next = mindex-DX;
  if (x>0)
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Added %d (%d, %d) with value %f.\n", next, x-1, y, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Updated %d (%d,%d) from %f to %f.\n",next, x-1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else					/* status is known. */
     {
      bpsi = psi[next];
      blabel = labels[next];
     }

    if (status[next] == KNOWN)
      lab_lf = labels[next];

    if (y>0)
     {
      next-=DY;
      if (status[next] == KNOWN)
       {
        lab_ul = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
    if (y<(NROWS-1))
     {
      next+=(DY+DY);
      if (status[next] == KNOWN)
       {
        lab_dl = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
   }

  next = mindex+DX;
  if (x<(NCOLS-1))
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Added %d (%d, %d) with value %f.\n",next, x+1, y, T);
      #endif

     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Updated %d (%d,%d) from %f to %f.\n",next, x+1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN) 
      lab_rt = labels[next];

    if (y<(NROWS-1))
     {
      next+=DY;
      if (status[next] == KNOWN)
       {
        lab_dr = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
    if (y>0)
     {
      next-=(DY+DY);
      if (status[next] == KNOWN)
       {
        lab_ur = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
   }

  next = mindex-DY;
  if (y>0)
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",next, x, y-1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Updated %d (%d,%d) from %f to %f.\n",next, x, y-1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN)
      lab_up = labels[next];
   }

  next = mindex+DY;
  if (y<(NROWS-1))
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", next, x, y+1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, NUMBER(1));

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Updated %d (%d,%d) from %f to %f.\n",next, x, y+1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN)
      lab_dn = labels[next];
   }

  if ( (status[mindex] != KNOWN) && blabel )
   {
    labels[mindex] = blabel;

    if ( lab_lf && (blabel != lab_lf) )
      shock[mindex] = 1;
    else if ( lab_rt && (blabel != lab_rt) )
      shock[mindex] = 2;
    else if ( lab_up && (blabel != lab_up) )
      shock[mindex] = 3;
    else if ( lab_dn && (blabel != lab_dn) )
      shock[mindex] = 4;
    else
      shock[mindex] = 0;
    // Using the extra if statements below will thicken the shock line.
    /*else if ( lab_ul && (blabel != lab_ul) )
      shock[mindex] = 5;
    else if ( lab_dl && (blabel != lab_dl) )
      shock[mindex] = 6;
    else if ( lab_ur && (blabel != lab_ur) )
      shock[mindex] = 7;
    else if ( lab_dr && (blabel != lab_dr) )
      shock[mindex] = 8;*/
   }
  else if (status[mindex] != KNOWN)
    labels[mindex] = 0;

  status[mindex] = KNOWN;
 } 

}

void fastmarch::labelmarch(number maxdist)
{
register number T;
register int next;
int mindex, x, y;
unsigned char testH, testV, testD1, testD2;
heapType el;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

//int iloop = minheap.numel();
//while (iloop-- > 0) 
while (minheap.pull(el))
 { 
  mindex = el.getIndex();

  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed (%d, %d) value %f.\n", ++cnt, x, y, 
           el.getTime());
  #endif


  if (psi[mindex] < maxdist)
   {
    x = XSUB(mindex);
    y = YSUB(mindex);
  
    testH = 0x00;
    testV = 0x00;
    testD1 = 0x00;
    testD2 = 0x00;
  
  
    next = mindex-DX;
    if (x>0)
     {
      if (status[next] == FAR)			/* add to heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
        status[next] = TRIAL;
        psi[next] = T;
        el.set(T,next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 2)
          DEBUGF("\t1--Added %d (%d, %d) with value %f.\n", next, x-1, y, T);
        #endif
       }
      else if (status[next] == TRIAL)		/* update heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t1--Updated %d (%d,%d) from %f to %f.\n",next, x-1, y, 
  	       psi[next], T);
        #endif
  
        psi[next] = T;
        minheap.update(minheap.getId(next),T);
       }
  
      if (status[mindex] != KNOWN)
       {
        if (psi[next] == psi[mindex])
          testH = 1;
        else if (psi[next] < psi[mindex])
          testH = 2;
       }
  
      if (y>0)
       {
        next-=DY;
        if (psi[next] == psi[mindex])
          testD1 = 1;
        else if (psi[next] < psi[mindex])
          testD1 = 2;
       }
      if (y<(NROWS-1))
       {
        next+=(DY+DY);
        if (psi[next] == psi[mindex])
          testD2 = 1;
        else if (psi[next] < psi[mindex])
          testD2 = 2;
       }
     }
  
  
    next = mindex+DX;
    if (x<(NCOLS-1))
     {
      if (status[next] == FAR)			/* add to heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
        status[next] = TRIAL;
        psi[next] = T;
        el.set(T,next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 2)
          DEBUGF("\t2--Added %d (%d, %d) with value %f.\n",next, x+1, y, T);
        #endif
  
       }
      else if (status[next] == TRIAL)		/* update heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t2--Updated %d (%d,%d) from %f to %f.\n",next, x+1, y, 
  	       psi[next], T);
        #endif
  
        psi[next] = T;
        minheap.update(minheap.getId(next),T);
       }
  
      if (testH)
       {
        if ( (testH == 1) && (psi[next] < psi[mindex]) )
          testH = 5;
        else if ( (testH == 2) && (psi[next] <= psi[mindex]) )
          testH = 5;
       }
  
      if (y<(NROWS-1))
       {
        next+=DY;
        if ( (testD1 == 1) && (psi[next] < psi[mindex]) )
          testD1 = 5;
        else if ( (testD1 == 2) && (psi[next] <= psi[mindex]) )
          testD1 = 5;
       }
      if (y>0)
       {
        next-=(DY+DY);
        if ( (testD2 == 1) && (psi[next] == psi[mindex]) )
          testD2 = 5;
        else if ( (testD2 == 2) && (psi[next] < psi[mindex]) )
          testD2 = 5;
       }
     }
  
    next = mindex-DY;
    if (y>0)
     {
      if (status[next] == FAR)			/* add to heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
        status[next] = TRIAL;
        psi[next] = T;
        el.set(T,next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 2)
          DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",next, x, y-1, T);
        #endif
       }
      else if (status[next] == TRIAL)		/* update heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t3--Updated %d (%d,%d) from %f to %f.\n",next, x, y-1,
                 psi[next], T);
        #endif
  
        psi[next] = T;
        minheap.update(minheap.getId(next),T);
       }
  
      if (status[mindex] != KNOWN)
       {
        if (psi[next] == psi[mindex])
          testV = 1;
        else if (psi[next] < psi[mindex])
          testV = 2;
       }
     }
  
    next = mindex+DY;
    if (y<(NROWS-1))
     {
      if (status[next] == FAR)			/* add to heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
        status[next] = TRIAL;
        psi[next] = T;
        el.set(T,next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 2)
          DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", next, x, y+1, T);
        #endif
       }
      else if (status[next] == TRIAL)		/* update heap. */
       {
        T = SOLVEDIST(next, NUMBER(1));
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t4--Updated %d (%d,%d) from %f to %f.\n",next, x, y+1,
                 psi[next], T);
        #endif
  
        psi[next] = T;
        minheap.update(minheap.getId(next),T);
       }
  
      if (testV)
       {
        if ( (testV == 1) && (psi[next] < psi[mindex]) )
          testV = 5;
        else if ( (testV == 2) && (psi[next] <= psi[mindex]) )
          testV = 5;
       }
     }
  
    shock[mindex] = (testH == 5) + ((testV == 5)<<1) 
                      + ((testD1 ==5)<<2) + ((testD2==5)<<3);
   } 
  status[mindex] = KNOWN;
 }

}
#undef LDEBUG

#undef NROWS
#undef NCOLS

   /*-------------------------- labelmarch -------------------------*/

/*-------------------------- labelfieldmarch -------------------------*/
/*



   Fast march out the label terms to identify the regions associated to 
   the different "sources." The traversal speed is given by a scalar field 
   (speed varies spatially).  While the heap is not empty, there are points 
   whose corresponding distances need to be updated.  Once the heap is empty, 
   there is no more work to do and the traversal time map should be bonafide.

   Inputs: N/A.

   Outputs: N/A.

   Requirements:
     -the traversal speed map should have been prepared initially be initDist.

   Modifies:

   Memory Usage: N/A.

   TODO: Does not generate a shock map like distmarch does.  Needs coding.
*/
#define NROWS M
#define NCOLS N

#define LDEBUG  0	/*  Local debug value.  Adjust to needs. */

void fastmarch::labelfieldmarch(number *speed)
{
register number T;
register int next;
int mindex, x, y;
labeltype lab_lf, lab_rt, lab_up, lab_dn;
labeltype lab_ul, lab_dl, lab_ur, lab_dr;
heapType el;

number bpsi;
labeltype blabel;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

while (minheap.pull(el))
 { 
  mindex = el.getIndex();
  x = XSUB(mindex);
  y = YSUB(mindex);

  blabel = 0;
  lab_lf = 0;
  lab_rt = 0;
  lab_up = 0;
  lab_dn = 0;
  lab_ul = 0;
  lab_dl = 0;
  lab_ur = 0;
  lab_dr = 0;

  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed (%d, %d) value %f.\n", ++cnt, x, y, 
           el.getTime());
  #endif

  next = mindex-DX;
  if (x>0)
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Added %d (%d, %d) with value %f.\n", next, x-1, y, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Updated %d (%d,%d) from %f to %f.\n",next, x-1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else					/* status is known. */
     {
      bpsi = psi[next];
      blabel = labels[next];
     }

    if (status[next] == KNOWN)
      lab_lf = labels[next];

    if (y>0)
     {
      next-=DY;
      if (status[next] == KNOWN)
       {
        lab_ul = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
    if (y<(NROWS-1))
     {
      next+=(DY+DY);
      if (status[next] == KNOWN)
       {
        lab_dl = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
   }

  next = mindex+DX;
  if (x<(NCOLS-1))
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Added %d (%d, %d) with value %f.\n",next, x+1, y, T);
      #endif

     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Updated %d (%d,%d) from %f to %f.\n",next, x+1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN) 
      lab_rt = labels[next];

    if (y<(NROWS-1))
     {
      next+=DY;
      if (status[next] == KNOWN)
       {
        lab_dr = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
    if (y>0)
     {
      next-=(DY+DY);
      if (status[next] == KNOWN)
       {
        lab_ur = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
   }

  next = mindex-DY;
  if (y>0)
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",next, x, y-1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Updated %d (%d,%d) from %f to %f.\n",next, x, y-1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN)
      lab_up = labels[next];
   }

  next = mindex+DY;
  if (y<(NROWS-1))
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", next, x, y+1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = SOLVEDIST(next, speed[mindex]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Updated %d (%d,%d) from %f to %f.\n",next, x, y+1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN)
      lab_dn = labels[next];
   }

  if ( (status[mindex] != KNOWN) && blabel )
   {
    labels[mindex] = blabel;

    if ( lab_lf && (blabel != lab_lf) )
      shock[mindex] = 1;
    else if ( lab_rt && (blabel != lab_rt) )
      shock[mindex] = 2;
    else if ( lab_up && (blabel != lab_up) )
      shock[mindex] = 3;
    else if ( lab_dn && (blabel != lab_dn) )
      shock[mindex] = 4;
    else
      shock[mindex] = 0;
    // Using the extra if statements below will thicken the shock line.
    /*else if ( lab_ul && (blabel != lab_ul) )
      shock[mindex] = 5;
    else if ( lab_dl && (blabel != lab_dl) )
      shock[mindex] = 6;
    else if ( lab_ur && (blabel != lab_ur) )
      shock[mindex] = 7;
    else if ( lab_dr && (blabel != lab_dr) )
      shock[mindex] = 8;*/
   }
  else if (status[mindex] != KNOWN)
    labels[mindex] = 0;

  status[mindex] = KNOWN;
 } 

}

#ifdef LDEBUG
#undef LDEBUG
#endif

#undef NROWS
#undef NCOLS

   /*----------------------- labelfieldmarch ----------------------*/

/*-------------------------- labelspeedmarch -------------------------*/
/*



   Fast march out the label terms to identify the regions associated to 
   the different "sources." The traversal speed is given by a scalar field 
   (speed varies spatially).  While the heap is not empty, there are points 
   whose corresponding distances need to be updated.  Once the heap is empty, 
   there is no more work to do and the traversal time map should be bonafide.

   Inputs: N/A.

   Outputs: N/A.

   Requirements:
     -the traversal speed map should have been prepared initially be initDist.

   Modifies:

   Memory Usage: N/A.

   TODO: Does not generate a shock map like distmarch does.  Needs coding.
*/
#define NROWS M
#define NCOLS N

#define LDEBUG  4	/*  Local debug value.  Adjust to needs. */

void fastmarch::labelspeedmarch(number *speed)
{
register number T;
register int next;
int mindex, x, y;
labeltype lab_lf, lab_rt, lab_up, lab_dn;
labeltype lab_ul, lab_dl, lab_ur, lab_dr;
labeltype lab_cn;
heapType el;

number bpsi;
labeltype blabel;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

while ((minheap.pull(el)) && (cnt < 1200))
 { 
  mindex = el.getIndex();
  x = XSUB(mindex);
  y = YSUB(mindex);
  lab_cn = labels[mindex]-1;

  blabel = 0;
  lab_lf = 0;
  lab_rt = 0;
  lab_up = 0;
  lab_dn = 0;
  lab_ul = 0;
  lab_dl = 0;
  lab_ur = 0;
  lab_dr = 0;

  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed %d - (%d, %d) value %f / label %d.\n", 
	       ++cnt, mindex, x, y, el.getTime(), lab_cn);
  #endif

  next = mindex-DX;
  if (x>0)
   {
    if (status[next] == FAR)				/* add to heap. */
     {
      T = SOLVEDIST(next, speed[lab_cn]);
	  labels[next] = labels[mindex];
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Added %d (%d, %d) with value %f.\n", next, x-1, y, T);
      #endif
     }
    else if ((status[next] == TRIAL) && (labels[next] == labels[mindex]))		
     { 										/* update heap. */
      T = SOLVEDIST(next, speed[lab_cn]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Updated %d (%d,%d) from %f to %f [%d, %d - %d].\n",
		       next, x-1, y, psi[next], T, labels[mindex], labels[next],
			   psi[next] > T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else									/* status is known. */
     {
      bpsi = psi[next];
      blabel = labels[next];
     }

    if (status[next] == KNOWN)
      lab_lf = labels[next];

    if (y>0)
     {
      next-=DY;
      if (status[next] == KNOWN)
       {
        lab_ul = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
    if (y<(NROWS-1))
     {
      next+=(DY+DY);
      if (status[next] == KNOWN)
       {
        lab_dl = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
   }

  next = mindex+DX;
  if (x<(NCOLS-1))
   {
    if (status[next] == FAR)				/* add to heap. */
     {
      T = SOLVEDIST(next, speed[lab_cn]);
	  labels[next] = labels[mindex];
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Added %d (%d, %d) with value %f.\n",next, x+1, y, T);
      #endif

     }
    else if ((status[next] == TRIAL) && (labels[next] == labels[mindex]))		
     { 										/* update heap. */
      T = SOLVEDIST(next, speed[lab_cn]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Updated %d (%d,%d) from %f to %f [%d, %d - %d].\n",
		       next, x+1, y, psi[next], T, labels[mindex], labels[next],
			   psi[next] > T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN) 
      lab_rt = labels[next];

    if (y<(NROWS-1))
     {
      next+=DY;
      if (status[next] == KNOWN)
       {
        lab_dr = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
    if (y>0)
     {
      next-=(DY+DY);
      if (status[next] == KNOWN)
       {
        lab_ur = labels[next];
        if ( !blabel || (bpsi > psi[next]) )
         {
          bpsi = psi[next];
          blabel = labels[next];
         }
       }
     }
   }

  next = mindex-DY;
  if (y>0)
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, speed[lab_cn]);
	  labels[next] = labels[mindex];
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",next, x, y-1, T);
      #endif
     }
    else if ((status[next] == TRIAL) && (labels[next] == labels[mindex]))		
     { 										/* update heap. */
      T = SOLVEDIST(next, speed[lab_cn]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Updated %d (%d,%d) from %f to %f [%d, %d - %d].\n",
		       next, x, y-1, psi[next], T, labels[mindex], labels[next],
			   psi[next] > T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN)
      lab_up = labels[next];
   }

  next = mindex+DY;
  if (y<(NROWS-1))
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      T = SOLVEDIST(next, speed[lab_cn]);
	  labels[next] = labels[mindex];
      status[next] = TRIAL;
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", next, x, y+1, T);
      #endif
     }
    else if ((status[next] == TRIAL) && (labels[next] == labels[mindex]))		
     { 										/* update heap. */
      T = SOLVEDIST(next, speed[lab_cn]);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Updated %d (%d,%d) from %f to %f [%d, %d - %d].\n",
		       next, x, y+1, psi[next], T, labels[mindex], labels[next],
			   psi[next] > T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
    else
     {
      if ( !blabel || (bpsi > psi[next]) )
       {
        bpsi = psi[next];
        blabel = labels[next];
       }
     }

    if (status[next] == KNOWN)
      lab_dn = labels[next];
   }

  if ( (status[mindex] != KNOWN) && blabel )
   {
    labels[mindex] = blabel;

    if ( lab_lf && (blabel != lab_lf) )
      shock[mindex] = 1;
    else if ( lab_rt && (blabel != lab_rt) )
      shock[mindex] = 2;
    else if ( lab_up && (blabel != lab_up) )
      shock[mindex] = 3;
    else if ( lab_dn && (blabel != lab_dn) )
      shock[mindex] = 4;
    else
      shock[mindex] = 0;
    // Using the extra if statements below will thicken the shock line.
    /*else if ( lab_ul && (blabel != lab_ul) )
      shock[mindex] = 5;
    else if ( lab_dl && (blabel != lab_dl) )
      shock[mindex] = 6;
    else if ( lab_ur && (blabel != lab_ur) )
      shock[mindex] = 7;
    else if ( lab_dr && (blabel != lab_dr) )
      shock[mindex] = 8;*/
   }
  else if (status[mindex] != KNOWN)
   {
	printf("Uh-Oh!!\n");
    labels[mindex] = 0;
   }

  status[mindex] = KNOWN;
 } 

}

#ifdef LDEBUG
#undef LDEBUG
#endif

#undef NROWS
#undef NCOLS

   /*----------------------- labelspeedmarch ----------------------*/

/*------------------------------ relabel -----------------------------*/
/*
   Perform an update of the extension data function to conform with 
   the new location of the zero-contour.  There are a variety of 
   extension data techniques available; the implementation chosen
   at compile time is selected via #define.

   Inputs: N/A.

   Outputs: N/A.

   Requirements: N/A.

   Modifies:
     -the functions, xphi, xphix, xphiy, will represent the proper extension
      data for the level set function.

   Memory Usage: 
     The _SIMPLE_RI_ and _PDE_RI_ implementations allocate a block of 
     integers of size area.  The block is de-allocated prior to function 
     termination.
*/

#define LDEBUG 1

#define NROWS M
#define NCOLS N

void fastmarch::relabel(void)
{

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::relabel -- missingInit",
         "Need to first initialize/define the initial label map.\n");
  return;
 }

compZeroContour(); 
initDist();
labelmarch();

}

#undef LDEBUG

#undef NROWS
#undef NCOLS

   /*--------------------------- relabel --------------------------*/


/*------------------------------ initNMS -----------------------------*/
/*
*/
void fastmarch::initNMS(number *values)
{
register int i, index;
int x, y, next;
heapType el;

/*------------- Initialize data. -------------*/

minheap.empty();				/* Prepare the heap. */

if (seedlen)					/* Done only if zero contour exists. */
 {
  for(i=area-1; i>=0; i--)		/* Clean slate, all pts FAR at INFINITY, */
   {							/*   and no shocks.                      */
    status[i] = FAR;
    psi[i] = 0;					/* psi encodes for the value of the max. */
    shock[i] = 0;
   }

  for(i=seedlen-1; i>=0; i--)	/* Start to fill in slate w/known data.  */
   {
    index = seed[i];
    status[index] = KNOWN;
	psi[index] = values[index];

    el.set(-values[index],index);	/* Place seed point into heap. */
    minheap.add(el);
   }
 }
else 
 {
  if (!ISNULL(labels))
   {
    for (i=area-1; i>=0; i--)
	 {
	  if (labels[i] > 0)
	   {
        status[i] = KNOWN;
		psi[i] = values[i];

        el.set(-values[i],i);	/* Place seed point into heap. */
        minheap.add(el);
	   }
	  else
	   {
        status[i] = FAR; 
        psi[i] = 0;
	   }

      shock[i] = 0;
     }
   }
 }

}
   /*--------------------------- initNMS --------------------------*/


/*-------------------------- compRegionsNMS --------------------------*/
/*
*/

void fastmarch::compRegionsNMS(number *values, number *threshold)
{
int i;

if ( ISNULL(labels) )
 {
  ERRORF("fastmarch::relabel -- missingInit",
         "Need to first initialize/define the initial label map.\n");
  return;
 }

initNMS(values);
nmsmarch(values, threshold);
}  
   /*---------------------- compRegionsNMS ---------------------*/

/*----------------------------- nmsmarch -----------------------------*/
/*

   March out the initial seed points descending in cost/value until
   it is no longer possible to do so.  Traversal speed is uniform, just
   expand out based on neighbors.

   Inputs: 
     

   Outputs: N/A.

   Requirements:
     - the seed points should be set.

   Modifies:

   Memory Usage: N/A.
     - 

*/
#define NROWS M
#define NCOLS N

#define LDEBUG  0	/*  Local debug value.  Adjust to needs. */

void fastmarch::nmsmarch(number *value, number *thresholds)
{
register int next;
int mindex, x, y;
heapType el;
number thresh;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

while (minheap.pull(el))
 { 
  mindex = el.getIndex();
  x = XSUB(mindex);
  y = YSUB(mindex);

  thresh = thresholds[labels[mindex]-1];

  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed (%d, %d) value %f.\n", ++cnt, x, y, 
           el.getTime());
  #endif

  next = mindex-DX;
  if (x>0)
   {
    if (status[next] == FAR) 
     {
	  if (value[next] > 0) 					/* add to heap. */
	   {
        status[next] = TRIAL;
        if (value[next] >= thresh)
          labels[next] = labels[mindex];
		else 
		  labels[next] = 0;
		psi[next] = psi[mindex];			/* keep track of the peak. */
        el.set(-value[next],next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t1a--Added %d (%d, %d) with value %f.\n", 
		         next, x-1, y, -value[next]);
        #endif
       }
      else
       status[next] = KNOWN;
     }
	else if ( (status[next] == TRIAL) && (psi[mindex] > psi[next]) )
	 {
	  labels[next] = labels[mindex];
	  psi[next] = psi[mindex];
	 }
   }

  next = mindex+DX;
  if (x<(NCOLS-1))
   {
    if (status[next] == FAR) 
	 {
      if (value[next] > 0)						/* add to heap. */
       {
        status[next] = TRIAL;
        if (value[next] >= thresh)
          labels[next] = labels[mindex];
		else 
		  labels[next] = 0;
		psi[next] = psi[mindex];			/* keep track of the peak. */
        el.set(-value[next],next);
        minheap.add(el);
  
        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t2a--Added %d (%d, %d) with value %f.\n",
		         next, x+1, y, -value[next]);
        #endif
       }
      else
       status[next] = KNOWN;
     }
	else if ( (status[next] == TRIAL) && (psi[mindex] > psi[next]) )
	 {
	  labels[next] = labels[mindex];
	  psi[next] = psi[mindex];
	 }
   }

  next = mindex-DY;
  if (y>0)
   {
    if (status[next] == FAR) 
	 {
      if (value[next] > 0)						/* add to heap. */
       {
        status[next] = TRIAL;
        if (value[next] >= thresh)
          labels[next] = labels[mindex];
		else 
		  labels[next] = 0;
		psi[next] = psi[mindex];			/* keep track of the peak. */
        el.set(-value[next],next);
        minheap.add(el);

        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",
		         next, x, y-1, -value[next]);
        #endif
       }
      else
       status[next] = KNOWN;
     }
	else if ( (status[next] == TRIAL) && (psi[mindex] > psi[next]) )
	 {
	  labels[next] = labels[mindex];
	  psi[next] = psi[mindex];
	 }
   }

  next = mindex+DY;
  if (y<(NROWS-1))
   {
    if (status[next] == FAR) 
     {
	  if (value[next] > 0)						/* add to heap. */
	   {
        status[next] = TRIAL;
        if (value[next] >= thresh)
          labels[next] = labels[mindex];
		else 
		  labels[next] = 0;
		psi[next] = psi[mindex];			/* keep track of the peak. */
        el.set(-value[next],next);
        minheap.add(el);

        #ifdef _DEBUG_
        if (LDEBUG > 3)
          DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", 
		         next, x, y+1, -value[next]);
        #endif
	   }
	  else
	   status[next] = KNOWN;
     }
	else if ( (status[next] == TRIAL) && (psi[mindex] > psi[next]) )
	 {
	  labels[next] = labels[mindex];
	  psi[next] = psi[mindex];
	 }
   }

  status[mindex] = KNOWN;
 } 

}

#ifdef LDEBUG
#undef LDEBUG
#endif

#undef NROWS
#undef NCOLS

   /*-------------------------- nmsmarch --------------------------*/

/* } */
/*========================== Label Propogation =========================*/


/*========================== Data Propogation ==========================*/
/* { */

/*------------------------------ setData -----------------------------*/
/*
  Sets the data of the fast march class to be the data argument.
  The data is copied, since the fast march class will extend the
  the data (unless, the _NOCOPY_ flag is set).  
  
  Input:
    image       - a pointer to the data.
    iM, iN, iP	- data dimensions.

  Output:
    N/A

  Changes:
    	      	- image points to the same memory location as 
                  img does.
*/

void fastmarch::setData(datatype *thedata, int iM, int iN, int iP)
{

if (!checkDimensions(iM, iN)) 
 {
  P = iP;
  volume = area * P;

  #ifdef _NOCOPY_
  xdata = thedata;
  #else
  if (ISNULL(thedata) && NOTNULL(xdata))
   {
    ::FREE(xdata);
    xdata = (datatype*)NULL;
   }
  else if (NOTNULL(thedata))
   {
    DEBUGF("Allocating data for extensions.\n");
    xdata = (datatype*)MALLOC(volume*sizeof(datatype));

    if (ISNULL(xdata))
     {
      ERRORF("fastmarch::setData -- allocProb",
             "Can't allocate memory for extended data.\n");
      return;
     }
    memcpy(xdata, thedata, volume*sizeof(datatype));
   }
  #endif
 }
else
 {
  ERRORF("fastmarch::setImage -- dim Error",
         "Dimensions of data incompatible with expected dimensions.");
  return;
 }

}  /*--------------------------- setData --------------------------*/


/*------------------------------ getData -----------------------------*/
/*
  Gets the data field from the fast march class.
  
  Input:
    image       - a pointer to the data.
    iM, iN, iP	- data dimensions.

  Output:
    N/A

  Changes:
    	      	- image points to the same memory location as 
                  img does.
*/

void fastmarch::getData(datatype *thedata)
{
if (NOTNULL(thedata) && NOTNULL(xdata))
 {
  DEBUGF("Copying extended data!.\n");
  memcpy(thedata, xdata, volume*sizeof(datatype));
 }
}

/*------------------------- getDataDimensions ------------------------*/
/*
*/

void fastmarch::getDataDimensions(int &iM, int &iN, int &iP)
{

iM = M;
iN = N;
iP = P;

}  /*---------------------- getDataDimensions ---------------------*/


/*------------------------ checkDataDimensions -----------------------*/
/*
*/

bool fastmarch::checkDataDimensions(int sM, int sN, int sP)
{
return !( (sM == M)  && (sN == N) && (sP == P));
}  /*--------------------- checkDataDimensions --------------------*/


/*----------------------------- solveData ----------------------------*/
/* 
   Solve for T at index using given speed function. Returns the new
   time value for the fast marching Boundary Value problem.  This
   function is specific to the fast march problem which computes
   both distances and extension velocities.  Only 4-neighbors are 
   considered.

   Inputs:

   Outputs:

   Requirements:

   Modifies:

   Memory Usage: N/A.
*/

#define NROWS M
#define NCOLS N

number fastmarch::solveData(int index, number speed)
{
number T, Tup, Tdn, Tlf, Trt;		/* neighboring T values.  */
number Dup, Ddn, Dlf, Drt;		/* diffs in 4 directions. */
int srcH, srcV;				/* source direction?      */
number A, B, C;				/* to solve quad. eqn.    */
float Disc;				/* discriminant 	  */

int x, y;
int i, offset;

srcH = 0;				/* initially set to no    */
srcV = 0; 				/*   contributions.       */

T = psi[index];

x = XSUB(index);
y = YSUB(index);

if (x>0)
  Tlf = psi[index-DX];
else
  Tlf = psi[index];

if (x<(NCOLS-1))
  Trt = psi[index+DX];
else
  Trt = psi[index];

if (y>0)
  Tup = psi[index-DY];
else
  Tup = psi[index];

if (y<(NROWS-1))
  Tdn = psi[index+DY];
else
  Tdn = psi[index];

Dup = T-Tup; 				/* discrete derivatives.         */
Ddn = Tdn-T;
Dlf = T-Tlf;
Drt = Trt-T;

A = NUMBER(0); 				/* setup quadratic equations.    */
B = NUMBER(0);
C = -INV(TIMES(speed,speed));

if ((Dlf>0) || (Drt<0))			/* check diff's for upwind soln. */
 {
  A++;
  if (Dlf > -Drt)			/* use left contribution.  */
   {
    B -= TIM2(Tlf);
    C += TIMES(Tlf,Tlf);
    srcH = -DX;
   }
  else					/* use right contribution. */
   {
    B -= TIM2(Trt);
    C += TIMES(Trt,Trt);
    srcH = DX;
   }
 }

if ((Dup>0) || (Ddn < 0))
 {
  A++;
  if (Dup > -Ddn) 			/* use up contribution.   */
   {
    B -= TIM2(Tup);
    C += TIMES(Tup,Tup);
    srcV = -DY; 		
   }
  else 					/* use down contribution. */
   {
    B -= TIM2(Tdn);
    C += TIMES(Tdn,Tdn);
    srcV = DY;
   }
 }

/* If all four differences were zero, A will be zero, and we return
   the original T value to avoid division by zero when solving the 
   quadratic equation.
*/
if (A!=0)
 {
  Disc = TIMES(B,B) - TIMES(TIMES(A,C),4);
  if (Disc >= 0)
   {
    /* Solve for T and extension velocity (plus sign correct due to
       upwind nature of problem).  Extension velocity is average
       of neighbors velocities.
    */
    T = DIVIDE( -B+NUMBER(sqrt(Disc)) , TIM2(A) );
   }
  else
   {
    /* Solve for T and extension velocity (plus sign correct due to
       upwind nature of problem).  Extension velocity is average
       of neighbors velocities.
    */
    if (T > Tup)
      T = Tup;
    if (T > Tdn)
      T = Tdn;
    if (T > Trt)
      T = Trt;
    if (T > Tlf)
      T = Tlf;
    T = T+INV(speed);
   }

  offset = 0;
  for (i=P;i>0;i--,offset+=area)
   {
    xdata[index+offset] = 0;
    if (srcH != 0)
      xdata[index+offset] += xdata[index+offset+srcH];
    
    if (srcV != 0)
      xdata[index+offset] += xdata[index+offset+srcV];

    xdata[index+offset] = DIVIDE(xdata[index+offset], A);
    //TODO: Should probably be a weighted average based on speed of arrival.
   }
 }

return T;
}

#undef NROWS
#undef NCOLS

   /*-------------------------- solveData -------------------------*/


/*----------------------------- datamarch ----------------------------*/
/*
  This function fast marches out an initialized distance map, where the
  traversal speed is given by a scalar field (speed varies spatially).  
  While the heap is not empty, there are points whose corresponding 
  distances need to be updated.  Once the heap is empty, there is no 
  more work to do and the traversal time map should be bonafide.

  Inputs: N/A.
  Outputs: N/A.
  Requirements:
    -the traversal speed map should have been prepared initially be initDist.

  Memory Usage: N/A.

*/
#define NROWS M
#define NCOLS N

#define LDEBUG 0	/*  Local debug value.  Adjust to needs. */

void fastmarch::datamarch(void)
{
register number T;
register int next;
int mindex, x, y;
heapType el;

#ifdef _DEBUG_
int cnt;
cnt = 0;

if (LDEBUG > 1)
  minheap.print();
#endif

while (minheap.pull(el))
 { 
  mindex = el.getIndex();
  status[mindex] = KNOWN;
  x = XSUB(mindex);
  y = YSUB(mindex);

  #ifdef _DEBUG_
  if (LDEBUG > 2)
    DEBUGF("Iteration %d grabbed (%d, %d) value %f.\n", ++cnt, x, y, 
           el.getTime());
  #endif

  next = mindex-DX;
  if (x>0)
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = solveData(next, (number)1);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t1--Added %d (%d, %d) with value %f.\n", next, x-1, y, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = solveData(next, (number)1);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t1--Updated %d (%d,%d) from %f to %f.\n",next, x-1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }

  next = mindex+DX;
  if (x<(NCOLS-1))
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = solveData(next, (number)1);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t2--Added %d (%d, %d) with value %f.\n",next, x+1, y, T);
      #endif

     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = solveData(next, (number)1);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t2--Updated %d (%d,%d) from %f to %f.\n",next, x+1, y, 
	       psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }

  next = mindex-DY;
  if (y>0)
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = solveData(next, (number)1);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t3--Added %d (%d, %d) with value %f.\n",next, x, y-1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = solveData(next, (number)1);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t3--Updated %d (%d,%d) from %f to %f.\n",next, x, y-1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }

  next = mindex+DY;
  if (y<(NROWS-1))
   {
    if (status[next] == FAR)			/* add to heap. */
     {
      status[next] = TRIAL;
      T = solveData(next, (number)1);
      psi[next] = T;
      el.set(T,next);
      minheap.add(el);

      #ifdef _DEBUG_
      if (LDEBUG > 2)
        DEBUGF("\t4--Added %d (%d, %d) with value %f.\n", next, x, y+1, T);
      #endif
     }
    else if (status[next] == TRIAL)		/* update heap. */
     {
      T = solveData(next, (number)1);

      #ifdef _DEBUG_
      if (LDEBUG > 3)
        DEBUGF("\t4--Updated %d (%d,%d) from %f to %f.\n",next, x, y+1,
               psi[next], T);
      #endif

      psi[next] = T;
      minheap.update(minheap.getId(next),T);
     }
   }
 } 

}

#undef LDEBUG

#undef NROWS
#undef NCOLS

   /*-------------------------- datamarch -------------------------*/


/*---------------------------- extendData ----------------------------*/
/*
*/
void fastmarch::extendData(void)
{

if ( ISNULL(xdata) )
 {
  ERRORF("fastmarch::relabel -- missingInit",
         "Need to first initialize/define the initial label map.\n");
  return;
 }

compZeroContour(); 
initDist();
datamarch();

}  /*------------------------- extendData -------------------------*/


/* } */
/*========================== Data Propogation ==========================*/


/*---------------------------- redistance ----------------------------*/
/*
   Perform redistancing.  There are a variety of redistancing techniques
   available; the implementation compiled is selected via #define.

   Inputs: N/A.
   Outputs: N/A.

   Modifies:
     -the distance map, psi, has been reinitialized.

   Requirements:

   Memory Usage: 
     The _SIMPLE_RI_ implementation allocates a block of integers of size
     area.  The block is de-allocated prior to function termination.
*/

#define LDEBUG 1

void fastmarch::redistance(void)
{

#if defined(_SIMPLE_RI_)
/*-- Simple one, preserve layer in and out of zero contour, then
     fast march out the new distances with respect to the preserved
     layers.  If a pixel is actually on the zero contour, then just
     use the pixel and not the laeyer in and out method.
--*/

int *tzero, tlen;	/* temporary zero contour and length. */
int *ozero, olen;	/* old zero contour and length. */
int i, x, y, marked;

tzero = (int *)CALLOC(area,sizeof(int));	/* allocate space (. */
if (!tzero)
 {
  ERRORF("levelset::redistance -- allocFail",
         "Failed to allocate temporary storage.");
  return;
 }

compZeroRegion(tzero, tlen);

ozero = seed;		/* Swap tzero and zero. */
olen  = seedlen;	/* Swap lengths also.   */
seed = tzero;
seedlen = tlen;

initDist();		/* Reinitialized distances by fast marching. */
distmarch();

seed = ozero;		/* Reset to original values. */
seedlen = olen;

::FREE(tzero);					/* deallocate space ). */

#elif defined(_INTERP_RI_)
/*-- Interpolate distance to front using piece-wise linear 
     approximation of the level set points that are near
     the boundary.  Points nears the boundary are those with
     neighbors who are on the other side of the zero contour
     (e.g., have different sign.).
--*/

interpDist();
distmarch();

#elif defined(_PDE_RI_)
/*-- Complex one, run PDE to generate new distances.  --*/

#ifdef _PDE_RI_HO_
distanceUpdatePDEHigherOrder();
#else
distanceUpdatePDE();
#endif

#else
/*-- Naive re-initialization of levelset.  Find pixel location of
     zero contour, then re-initialize level set.   This actually
     sets the "crudely" approximated zero contour distances to be
     zero, then fast marches out the distances relative to this
     contour.
--*/

//TODO: Need to correct.
compZeroContour(); 		
initDist();
distmarch();


#endif

}

#undef LDEBUG

   /*------------------------- redistance -------------------------*/



/*-------------------------- forceRedistance -------------------------*/
/*
*/

void fastmarch::forceRedistance(void)
{
redistance();
}  /*----------------------- forceRedistance ----------------------*/


/*------------------------- compdistCostField ------------------------*/
/*
    Memory:
      Allocates space to convert desired cost field into an effective
      cost field for the fast marching function.  The presumed
      relationship is inverse.  If other monotonic relationships are
      desired for conversion between cost and speed, then conversion should
      be done beforehand with compdistSpeedField as the proper marching
      function.
*/

void fastmarch::compdistCostField( number *cost )
{
int i;
number *speed;

speed = (number *)MALLOC(area*sizeof(number));
if (ISNULL(speed))
 {
  ERRORF("levelset::compdistCostField -- allocFail",
         "Could not allocate return space.");
 }
else
 {
  for (i=area-1;i>=0;i--)
    speed[i] = 1/cost[i];
 }

initDist();
fieldmarch(speed);
}  /*---------------------- compdistCostField ---------------------*/

/*------------------------ compdistSpeedField ------------------------*/
/*
*/

void fastmarch::compdistSpeedField( number *speed )
{
int i;

initDist();
fieldmarch(speed);
}  /*--------------------- compdistSpeedField ---------------------*/


/*--------------------------- getDistance ----------------------------*/
/*

   Memory Usage:
     Allocates return space for Distance function, psi.  Must be freed 
     manually outside of this function.
*/
number* fastmarch::getDistance(void)
{
number *retPsi;
register int i;

retPsi = (number *)MALLOC(area*sizeof(number));

if (!retPsi)
 {
  ERRORF("levelset::getPsi -- allocFail",
         "Could not allocate return space.");
 }
else
 {
  for (i=area-1;i>=0;i--)
    retPsi[i] = ((number)sign[i])*psi[i];
 }

return retPsi;
}

void fastmarch::getDistance(number *rpsi)
{
register int i;

for (i=area-1;i>=0;i--)
  rpsi[i] = ((number)sign[i])*psi[i];
}
   /*------------------------ getDistance -------------------------*/


/*--------------------------- getShockMap ----------------------------*/
/*

   Memory Usage:
     Allocates return space for shock map, rshock.  Must be freed 
     manually outside of this function.
*/
shocktype* fastmarch::getShockMap(void)
{
shocktype *rshock;
register int i;

rshock = (shocktype *)MALLOC(area*sizeof(shocktype));

if (ISNULL(rshock))
 {
  ERRORF("levelset::getPsi -- allocFail",
         "Could not allocate return space.");
 }
else
  memcpy(rshock, shock, sizeof(shocktype)*area);

return rshock;
}

void fastmarch::getShockMap(shocktype *rshock)
{

memcpy(rshock, shock, sizeof(shocktype)*area);

}
   /*------------------------ getShockMap -------------------------*/



/*------------------------------- free -------------------------------*/

#define LDEBUG 0

void fastmarch::free(void)
{

/*-------- De-allocate levelset memory. --------*/
if NOTNULL(psi)
 {
  ::FREE(psi);
  psi = (number*)NULL;
 }
if NOTNULL(sign)
 {
  ::FREE(sign);
  sign = (signtype*)NULL;
 }
if NOTNULL(phi)
 {
  ::FREE(phi);
  phi = (number *)NULL;
 }


if NOTNULL(labels)
 {
  ::FREE(labels);
  labels = (labeltype *)NULL;
 }
if NOTNULL(shock)
 {
  ::FREE(shock);
  shock = (shocktype *)NULL;
 }
#ifndef _NOCOPY_
if NOTNULL(xdata)
 { 
  ::FREE(xdata);
  xdata = (datatype*)NULL;
 }
#endif

if NOTNULL(status) 
 {
  ::FREE(status);
  status = (fmstat*)NULL;
 }
if NOTNULL(seed)
 {
  ::FREE(seed);
  seed = (int *)NULL;
 }
seedlen = 0;

minheap.free();

#ifdef _MEX_
if ( isPersistent && NOTNULL(matVar) )
 {
  ARRAYFREE(matVar);
  matVar = (mxArray*)NULL;
 }
#endif

#ifdef _DEBUG_
if ((DEBUG > 2) || (LDEBUG > 1))
  DEBUGF("\tDone de-allocating level set components.\n");
#endif

}  

#undef LDEBUG

   /*---------------------------- free ----------------------------*/



/*---------------------------- ~fastmarch ----------------------------*/

fastmarch::~fastmarch(void)
{
}  /*------------------------- ~fastmarch -------------------------*/



/*============================== fastmarch =============================*/

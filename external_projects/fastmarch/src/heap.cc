/*================================ heap ================================*/
/*
  Author:	Patricio Vela (pvela@ece.gatech.edu)
  Created:	02/02/2004
  Modified:	02/02/2004

  Implements a minHeap.  The principal property of a minHeap is
  that the root node represents the minimum element of all elements
  in the heap.  Furthermore, the children of any node are less than 
  or equal to the parent node.  The elements of a minHeap must be
  part of an ordered set with defined binary operations of <, >, and =.

  With a minHeap, the minimum element can be found in O(1) operations,
  the addition, removal, or modification of elements involves an
  update which takes O(log N) operations.

  History:

   ***** v2.3 - 08/12/2005 *****
   ***** v2.2 - 05/09/2005 *****
   ***** v2.1 - 03/25/2005 *****
   ***** v2.0 - 10/07/2004 *****
   ***** v1.1 - 03/24/2004 *****
   *
   * -modified to be specific to the levelset implementation.
   * -class no longer templated; code for heapType class added.
   * -attempting to integrate MATLAB mex support (not fully done).
   *
   *****

   ***** v1.0 - 02/02/2004 *****
   *
   * -source file created.
   * -heap class is templated.
   *
   *****

   Copyright(c) 2005.
*/
/*=============================== heap ===============================*/
#define __HEAP_CC

#include<number.h>
#include<heap.h>

#include<stdlib.h>

#ifdef _MEX_
# include<mex.h>
#else
# ifdef _DEBUG_
#  include<stdio.h>
# endif
#endif

#include<mymex.h>
#include<extras.h>



#define DEBUG 0

/*============================= heapType =============================*/


/*----------------------------- heapType -----------------------------*/
/*

*/

heapType::heapType(void)
{
}

/*---------------------------- ~heapType -----------------------------*/
/*

*/

heapType::~heapType(void)
{
}

/*-------------------------------- set -------------------------------*/
/*

*/

void heapType::set(number stime, int sindex)
{
time  = stime;
index = sindex;
}

/*------------------------------ setTime -----------------------------*/
/*

*/

void heapType::setTime(number ntime)
{
time = ntime;
}

/*------------------------------ getTime -----------------------------*/
/*

*/

number heapType::getTime(void)
{
return time;
}


/*----------------------------- setIndex -----------------------------*/
/*

*/
void heapType::setIndex(int nindex)
{
index = nindex;
}


/*----------------------------- getIndex -----------------------------*/
/*

*/
int heapType::getIndex(void)
{
return index;
}


/*---------------------------- operator<  ----------------------------*/
/*
  
*/

int heapType::operator<(heapType &other)
{
return (time < other.time);
}

int heapType::operator<(number ntime)
{
return (time < ntime);
}

/*---------------------------- operator>  ----------------------------*/
/*

*/

int heapType::operator>(heapType &other)
{
return (time > other.time);
}

int heapType::operator>(number ntime)
{
return (time > ntime);
}

/*---------------------------- operator==  ----------------------------*/
/*

*/

int heapType::operator==(heapType &other)
{
return (time == other.time);
}

int heapType::operator==(number ntime)
{
return (time == ntime);
}

/*---------------------------- operator=  ----------------------------*/
/*

*/

number heapType::operator=(heapType &other)
{
time  = other.time;
index = other.index;

return (time);
}


/*============================= heapType =============================*/


/*=============================== heap ===============================*/

/*------------------------------- heap -------------------------------*/
/*

*/


heap::heap(void)
{

data = NULL;
back = NULL;
last = -1;

}


/*------------------------------- ~heap ------------------------------*/
/*

*/


heap::~heap(void)
{

free();

data = NULL;
back = NULL;
size = 0;
last = -1;

}


/*------------------------------- init -------------------------------*/
/*

*/


void heap::init(int hsize)
{

if (hsize >= 0)
 {
  if (data == NULL)
   {
    data = new heapType[hsize];
    back = (backType*)CALLOC(hsize, sizeof(backType));
    size = hsize;
  
    if (data == NULL)
      ERRORF("heap::new -- allocFail",
             "Could not allocate space for heap.");
   }
  else
    DEBUGF("heap:new, heap appears to exist already.");
 }
else
 DEBUGF( "heap::new, heap size should be non-negative.");

}

/*------------------------------- free -------------------------------*/
/*
   Frees/de-allocates all space assigned for use by the heap class.
*/

void heap::free(void)
{
void *tmp;

if (size > 0)
 {
  if (data == NULL)
   {
    if (DEBUG > 1)
      DEBUGF("heap::free, heap does not appear to exist.");
   }
  else
   {
    delete [] data;
    ::FREE((void *) back);
   }
 }
#ifdef _DEBUG_
else if (DEBUG)
  DEBUGF("heap::free, heap size means heap occupies no space.\n"); 
#endif

data = NULL;
back = NULL;
size = 0;
last = -1;

}

/*------------------------------ isEmpty -----------------------------*/
/*

*/


inline heapId heap::isEmpty(void)
{

return (last+1);

}

/*------------------------------- empty ------------------------------*/
/*

*/


void heap::empty(void)
{

last = -1;

}

/*-------------------------------- add -------------------------------*/
/*

*/


int heap::add(heapType& element)
{

if (last < (size-1))
 {
  data[++last] = element;
  back[element.getIndex()] = last;
  promote(last);
  return(1);
 }
else
 {
  DEBUGF("heap::add, cannot add items.  heap is full.");
  return(0);
 }

}

/*------------------------------- pull -------------------------------*/
/*

*/


int heap::pull(heapType& boss)
{

boss = data[0];

if (last > 0)
 {
  back[data[last].getIndex()] = 0;
  data[0] = data[last--];
  sift(0);
 }
else if (last == 0)
  last--;
else
  return(0);

return(1);
}


/*------------------------------ update ------------------------------*/
/*

*/


void heap::update(heapId id, heapType &item)
{
heapType old;

if (item > data[id])
 {
  data[id] = item;
  sift(id);
 }
else if (item < data[id])
 {
  data[id] = item;
  promote(id);
 }

}

void heap::update(heapId id, number ntime)
{
heapType old;

if (data[id] < ntime)
 {
  data[id].setTime(ntime);
  sift(id);
 }
else if (data[id] > ntime)
 {
  data[id].setTime(ntime);
  promote(id);
 }

} /*--------------------------- update ---------------------------*/

/*------------------------------- sift -------------------------------*/
/*
   The sift routine looks to the bottom of the heap to ensure that the 
   indexed heap item is indeed smaller than those below it.  The heap 
   property is garaunteed only if the indexed heap item is smaller that 
   those above it.
*/

void heap::sift(heapId ind)
{
heapId left, right, minid;
heapType tmp;
int tback;

if (ind > 0)
 {
  left  = (ind<<1) + 1;
  right = (ind<<1) + 2;
 }
else
 {
  left  = 1;
  right = 2;
 }

while (left <= last)			// while not at bottom of heap:
 {
  minid = left;				// smallest child is left one.
  
  if (right <= last)			// if there's a child to right,
    if (data[right] < data[left])	// check to see if smaller,
      minid = right;			// to grab index as needed.

  if (data[minid] < data[ind])		// if child is smaller than leaf,
   {
    tmp = data[ind];
    tback = back[tmp.getIndex()];

    data[ind] = data[minid];		// then swap data, 
    back[tmp.getIndex()] = back[data[minid].getIndex()];

    back[data[minid].getIndex()] = tback;
    data[minid] = tmp;

    ind = minid;			// update current leaf to sift, 
    left  = (ind<<1) + 1;		// and
    right = (ind<<1) + 2;		// update location of children.
   }
  else
    break;				// otherwise, quit loop.

 }

}  /*---------------------------- sift ----------------------------*/

/*------------------------------ promote -----------------------------*/
/*

*/


void heap::promote(heapId ind)
{
heapId parent;
heapType tmp;
int tback;

while (ind > 0)				// while not at top of heap:
 {
  parent  = (ind-1)>>1;			// get location of parent node.
  
  if (data[ind] < data[parent])		// if child smaller than parent,
   {
    tmp = data[ind];
    tback = back[data[ind].getIndex()];

    data[ind] = data[parent];		// then swap data, and
    back[tmp.getIndex()] = back[data[parent].getIndex()];

    back[data[parent].getIndex()] = tback;
    data[parent] = tmp;

    ind = parent;			// update curr. leaf to promote. 
   }
  else
    break;				// otherwise, quit loop.
 }

}  /*-------------------------- promote --------------------------*/


/*------------------------------- getId ------------------------------*/
/*

*/


heapId heap::getId(int index)
{

return back[index];

}

#ifdef _MEX_

/*---------------------------- mexPersist ----------------------------*/
/*
*/

void heap::mexPersist()
{
if (NOTNULL(back))
  PERSIST(back);
} /*---------------------------- mexPersist ----------------------------*/
#endif // _MEX_


#ifdef _DEBUG_

/*------------------------------- numel ------------------------------*/
/*
*/

int heap::numel(void)
{
return  size;
}  /*---------------------------- numel ---------------------------*/

/*------------------------------- print ------------------------------*/
/*

*/


void heap::print(void)
{
heapId ind;

printf("\n");
for (ind=0; ind <= last; ind++)
#ifdef _INTEGER_
 printf("  index %d has heap val %d and index %d (back %d).\n",
        ind,data[ind].getTime(),data[ind].getIndex(),
	back[data[ind].getIndex()]);
#else
 printf("  index %d has heap val %5.3f and index %d (back %d).\n",
        ind,data[ind].getTime(),data[ind].getIndex(),
	back[data[ind].getIndex()]);
#endif

}  /*---------------------------- print ---------------------------*/

#endif /* _DEBUG_ */

/*================================ heap ================================*/

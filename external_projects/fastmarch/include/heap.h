/*=============================== heap ===============================*/
/*

  Author:	Patricio Vela (pvela@ece.gatech.edu)

  Created:	02/02/2004
  Modified:	08/13/2005

  Header file for a minHeap implementation.  
  This one is specific to a levelset implementation in that is 
  contains a back pointer to allow for quick indexing and modification 
  of updated values.  Furthermore, a class named heapType is defined.
  The elements of the heap are of type heapType.

  History:

   ***** v2.4 *****
   ***** v2.3 - 08/12/2005 *****
   ***** v2.2 - 05/09/2005 *****
   ***** v2.1 - 03/25/2005 *****
   ***** v2.0 - 10/07/2004 *****
   ***** v1.1 - 03/24/2004 *****
   *
   * -modified to be specific to levelset implementation.
   * -class definition no longer templated.
   *
   *****

   ***** v1.0 - 02/02/2004 *****
   *
   * -header file created.
   * -heap class is templated.
   *
   *****

   Copyright(c) 2005.
*/    
/*=============================== heap ===============================*/
#ifndef __HEAP_H
#define __HEAP_H


typedef int heapId;
typedef heapId backType;

/*-------------- Define heapType. --------------*/

class heapType
{

public:
  heapType(void);
  ~heapType(void);

  int operator<(heapType&);
  int operator>(heapType&);
  int operator==(heapType&);
  number operator=(heapType&);

  int operator<(number);
  int operator>(number);
  int operator==(number);

  void set(number stime,int sindex);
  void setTime(number ntime);
  number getTime(void);
  void setIndex(int nindex);
  int getIndex(void);

private:
  int index;
  number time;
};


/*------------- Define heap class. -------------*/

class heap
{

public:
  heap(void);			// constructor.
  ~heap(void);			// destructor.

  void init(int);			// initialize the heap.
  void free(void);			// free the heap.

  heapId isEmpty(void);			// is heap empty or not?
  void empty(void);			// empty the heap.
  int  add(heapType&);			// add element to heap.
  int  pull(heapType&);			// grab element from heap.
  void update(heapId, heapType&);	// update heap node data.
  void update(heapId, number);		// update heap node values.
  heapId getId(int);			// get Id from index.

  #ifdef _MEX_
    void mexPersist();
  #endif

  int numel(void);

  #ifdef _DEBUG_
    void print(void);
  #endif 

private:
  void sift(heapId);			// adjust heap if value increases.
  void promote(heapId);			// adjust heap if value decreases.

  heapType *data;	//  Pointer to the data (array of elements).
  backType *back;	//  Back-pointer, grid index -> heap location.
  heapId last;		//  Index to last element in heap.
  int size;		//  The total size of the heap.
};



#endif /* __HEAP_H */
/*=============================== heap ===============================*/

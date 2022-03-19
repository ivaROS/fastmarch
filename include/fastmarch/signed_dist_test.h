/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

/* distance transform */

#ifndef SIGNED_DT_H
#define SIGNED_DT_H

#include <algorithm>
#include <image.h>
#include <misc.h> 

#define INF 1E20

/* dt of 1d function using squared distance */
static float *dt(float *f, int n) {
  float *d = new float[n];
  int *v = new int[n];
  float *z = new float[n+1]; 
  int k = 0; //index of right-most parabola in the lower envelope
  v[0] = 0; //locations of parabolas in lower envelope
  z[0] = -INF; //locations of boundaries between parabolas. eg. [z[i],z[i+1]] is the range in which the i-th parabola of the lower envelope is below the others
  z[1] = +INF;
  for (int q = 1; q <= n-1; q++) {
    float s  = ((f[q]+square(q))-(f[v[k]]+square(v[k])))/(2*q-2*v[k]);  //horizontal position of the intersection between point at gridpoint q and right-most parabola
    while (s <= z[k]) { //if intersection is to the left of starting range of right-most parabola
      k--;  //remove right-most parabola
      s  = ((f[q]+square(q))-(f[v[k]]+square(v[k])))/(2*q-2*v[k]);
    }
    k++; //add new parabola
    v[k] = q; //add new parabola
    z[k] = s; //update previous parabola's range end
    z[k+1] = +INF; //set new parabola's range end
  }

  k = 0;
  for (int q = 0; q <= n-1; q++) {
    while (z[k+1] < q)  //set k to parabola whose range includes gridpoint q
      k++;
    d[q] = square(q-v[k]) + f[v[k]];  //compute squared distance for q
  }

  delete [] v;
  delete [] z;
  return d;
}


//row-wise signed distance
float *sdt(float *f, int n)
{
  //Analyze from left to right
  for(int q = 1; q < n; q++)
  {
    auto pv = f[q-1];
    auto& cv = f[q];
    
    if(pv>0)  //Previous was unoccupied
    {
      if(cv>0)  //Current is unoccupied
      {
        cv = pv+1;  //Increment distance
      }
      else
      {
        cv = 0; //Start of new occupied region
      }
    }
    else      //Previous was occupied
    {
      if(cv>0)  //Current is unoccupied
      {
        cv = 1;  //Reset distance
      }
      else
      {
        cv = pv-1;  //Increment negative distance
      }
    }
  }
  
  //Analyze from right to left
  for(int q = n-1; q >0; q--)
  {
    auto pv = f[q];
    auto& cv = f[q-1];
    
    if(pv>0)  //Previous was unoccupied
    {
      if(cv>0)  //Current is unoccupied
      {
        cv = std::min(cv, pv+1);  //the lesser of left and right distances
      }
      else
      {
        cv = 0; //Start of new occupied region
      }
    }
    else      //Previous was occupied
    {
      if(cv>0)  //Current is unoccupied
      {
        cv = 1;  //Reset distance
      }
      else
      {
        cv = std::max(cv, pv-1);  //the lesser of left and right (negative) distances
      }
    }
  }
  
  return f;
}



/* dt of 2d function using squared distance */
static void dt(image<float> *im) {
  int width = im->width();
  int height = im->height();
  float *f = new float[std::max(width,height)];

  // transform along columns
  for (int x = 0; x < width; x++) {
    //fill f with column entries
    for (int y = 0; y < height; y++) {
      f[y] = imRef(im, x, y);
    }
    //compute transform of f
    float *d = dt(f, height);
    //copy transformed results back into column
    for (int y = 0; y < height; y++) {
      imRef(im, x, y) = d[y];
    }
    delete [] d;
  }

  // transform along rows
  for (int y = 0; y < height; y++) {
    //fill f with row entries (why not just operate over image row directly?)
    for (int x = 0; x < width; x++) {
      f[x] = imRef(im, x, y);
    }
    //compute transform of f
    float *d = dt(f, width);
    //copy transformed results back into row
    for (int x = 0; x < width; x++) {
      imRef(im, x, y) = d[x];
    }
    delete [] d;
  }

  delete f;
}


/* dt of binary image using squared distance */
static image<float> *dt(image<uchar> *im, uchar on = 1) {
  int width = im->width();
  int height = im->height();

  //initialize distances of occupied points to 0 and the rest to a big number
  image<float> *out = new image<float>(width, height, false);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (imRef(im, x, y) == on)
        imRef(out, x, y) = -INF;
      else
        imRef(out, x, y) = INF;
    }
  }

  dt(out);
  return out;
}

#endif

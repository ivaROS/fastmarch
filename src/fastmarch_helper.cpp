#include <fastmarch/fastmarch_helper.h>

#include <opencv2/core/mat.hpp>
#include <opencv2/highgui.hpp>

/*------------- Includes -------------*/
#include<stdio.h>

#include<mymex.h>
#include<number.h>
#include<heap.h>
#include<extras.h>
#include<fastmarch.h>


/*------------ Definitions -----------*/

#define DEBUG 0

namespace fastmarching
{

/* TODO: Inherit from fastmarch (protected) in order to access lower level functions directly and avoid unnecessary recomputation.
 * For example: if labels have been computed w/ compLabels(), it is only necessary to run distmarch() to get distances
 */

  FastMarchHelper::FastMarchHelper():
    fmobj_(new fastmarch())
  {
    
  }
  
  FastMarchHelper::~FastMarchHelper()
  {
    fmobj_->free();
    delete fmobj_;
  }
  
  static cv::Mat convertMat(cv::Mat input, int type)
  {
    if(input.type() != type)
    {
      input.convertTo(input, type);
    }
    else if(!input.isContinuous())
    {
      input = input.clone();
    }
    return input;
  }
  
  
  cv::Mat FastMarchHelper::bwdist(cv::Mat binary_mask)
  {
    mM = binary_mask.rows;
    mN = binary_mask.cols;
    
    binary_mask = convertMat(binary_mask > 0, CV_64F);
    
    fmobj_->init(mM, mN);
    fmobj_->setSeedRegions((masktype*)binary_mask.data, mM, mN);
    fmobj_->forceRedistance();
    
    return getDistance();
  }
  
  cv::Mat FastMarchHelper::labelMarch(cv::Mat input_labels, cv::Mat* shockmap)
  {
    mM = input_labels.rows;
    mN = input_labels.cols;
    
    input_labels = convertMat(input_labels, CV_32S);
    
    fmobj_->init(mM, mN);
    fmobj_->setLabelRegions((labeltype *)input_labels.data, mM, mN);
    fmobj_->compLabels();
    
    if(shockmap)
    {
      *shockmap = getShockMap();
    }

    return getLabels();
  }
  
  
  cv::Mat FastMarchHelper::getLabels()
  {
    int type = CV_32S;
    cv::Mat result_labels(mM, mN, type);
    fmobj_->getLabels( (labeltype *)result_labels.data );
    return result_labels;
  }
  
  cv::Mat FastMarchHelper::getDistance()
  {
    int type = CV_32F;
    cv::Mat distance(mM, mN, type);
    fmobj_->getDistance( (number *)distance.data );
    return distance;
  }
  
  cv::Mat FastMarchHelper::getShockMap()
  {
    int type = CV_8U;
    cv::Mat shockmap(mM, mN, type);
    fmobj_->getShockMap( (shocktype *)shockmap.data );
    return shockmap;
  }
  
  
  void FastMarchHelper::showDistance()
  {
    cv::Mat distance = getDistance();
    cv::Mat distance_viz;
    cv::normalize(distance, distance_viz, 1, 0, cv::NORM_MINMAX);
    cv::imshow("distance", distance_viz);
  }
  
  void FastMarchHelper::showLabels()
  {
    cv::Mat result_labels = getLabels();
    cv::Mat result_labels_viz;
    result_labels.convertTo(result_labels_viz, CV_32FC1); // or CV_32F works (too)
    cv::normalize(result_labels_viz, result_labels_viz, 1, 0, cv::NORM_INF);
    cv::imshow("labels", result_labels_viz);
  }
  
} //end namespace fastmarch

#ifndef FASTMARCH_FASTMARCH_HELPER_H
#define FASTMARCH_FASTMARCH_HELPER_H

#include <opencv2/core/mat.hpp>

class fastmarch;

namespace fastmarching
{

  
  /* TODO: Inherit from fastmarch (protected) in order to access lower level functions directly and avoid unnecessary recomputation.
  * For example: if labels have been computed w/ compLabels(), it is only necessary to run distmarch() to get distances
  */
  class FastMarchHelper
  {
    fastmarch* fmobj_; 
    int mM, mN;
    
  public:
    FastMarchHelper();
    
    ~FastMarchHelper();
    
    cv::Mat bwdist(cv::Mat binary_mask);
    
    cv::Mat labelMarch(cv::Mat input_labels, cv::Mat* shockmap=nullptr);
    
    cv::Mat getLabels();
    
    cv::Mat getDistance();
    
    cv::Mat getShockMap();
    
    void showDistance();
    
    void showLabels();
    
  };

} //end namespace fastmarch

#endif //FASTMARCH_FASTMARCH_HELPER_H

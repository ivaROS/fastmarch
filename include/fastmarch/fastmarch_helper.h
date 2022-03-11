#ifndef FASTMARCH_FASTMARCH_HELPER_H
#define FASTMARCH_FASTMARCH_HELPER_H

#include <opencv2/core/mat.hpp>
#include <opencv2/highgui.hpp>  //The visualization capabilities should be moved to a separate class

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
    
    void init(cv::Mat mat);
    
    cv::Mat bwdist(cv::Mat binary_mask);
    
    cv::Mat labelMarch(cv::Mat input_labels, cv::Mat* shockmap=nullptr);
    
    cv::Mat getLabels() const;
    
    cv::Mat getDistance() const;
    
    cv::Mat getShockMap() const;
    
    void showDistance() const;
    
    void showLabels() const;
    
  };

} //end namespace fastmarch

#endif //FASTMARCH_FASTMARCH_HELPER_H

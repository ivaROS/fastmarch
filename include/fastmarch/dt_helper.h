#include <opencv2/core/mat.hpp>

namespace fastmarching
{
  
  class DTHelper
  {
    
  public:
    DTHelper();
    
    ~DTHelper();
    
    cv::Mat run(cv::Mat raw_mat);
    
  };

} //end namespace fastmarching

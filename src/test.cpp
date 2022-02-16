#include <fastmarch/fastmarch_helper.h>
#include <opencv2/highgui.hpp>


int main(void)
{
  fastmarching::FastMarchHelper fmh;
  cv::Mat labels = cv::Mat::zeros(100,100,CV_32S);
  labels.at<int>(50,50) = 1;
  labels.at<int>(50,51) = 2;

  fmh.bwdist(labels);
  fmh.showDistance();
  
  fmh.labelMarch(labels);
  fmh.showLabels();
  
  cv::waitKey(0);
  
  
  return 0;
}

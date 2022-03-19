#include <fastmarch/dt_helper.h>
#include <fastmarch/fastmarch_helper.h>
#include <opencv2/highgui.hpp>
#include <iostream>

int main(void)
{
  fastmarching::FastMarchHelper fmh;
  int dim=11;
  cv::Mat labels = cv::Mat::zeros(128,512,CV_32S);
  labels.at<int>(64,256) = 1;
  //labels.at<int>(dim/2,dim/2+1) = 2;

  fmh.bwdist(labels);
  fmh.showDistance();
  
  fmh.labelMarch(labels);
  fmh.showLabels();
  
  cv::Mat fmdist = fmh.getDistance();
  
  fastmarching::DTHelper dth;
  cv::Mat result = dth.run(labels);
  cv::Mat result_viz;
  cv::normalize(result, result_viz, 1, 0, cv::NORM_INF);
  cv::imshow("distance_transform", result_viz);
  
  cv::Mat diff = fmdist - result;
  
  cv::Mat diff_viz;
  cv::normalize(diff, diff_viz, 1, 0, cv::NORM_INF);
  cv::imshow("diff", diff_viz);
  
  std::cout << "Fastmarch:\n" << fmdist << "\n\n\nDT:\n" << result <<  "\n\n\nDiff:\n" << diff;

  cv::waitKey(0);
  
  return 0;
}

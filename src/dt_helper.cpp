#include <fastmarch/dt_helper.h>

#include <dt/misc.h>  //Should have been included by dt.h
#include <dt/dt.h>
#include <dt/imutil.h>
#include <dt/imconv.h>

#include <opencv2/core/mat.hpp>
#include <cstring>


namespace fastmarching
{

    DTHelper::DTHelper(){}
    
    DTHelper::~DTHelper(){}
    
    cv::Mat DTHelper::run(cv::Mat raw_mat)
    {
      cv::Mat mat = raw_mat>0;
      
      image<uchar> *input = new image<uchar>(mat.cols, mat.rows, false);
      for(int i = 0; i<mat.rows; i++)
      {
          const uchar* p = mat.ptr<uchar>(i);
          uchar* ip = input->access[i];
          for(int j = 0; j<mat.cols; j++)
          {
            ip[j] = p[j];
          }
      }
      
      image<float> *out = dt(input, 255);
      // take square roots
      for (int y = 0; y < out->height(); y++) {
        for (int x = 0; x < out->width(); x++) {
          imRef(out, x, y) = sqrt(imRef(out, x, y));
        }
      }
      
      cv::Mat out_mat(mat.cols, mat.rows, CV_32F);
      std::memcpy(out_mat.data, out->data, mat.cols*mat.rows*sizeof(float));
      
      delete input;
      delete out;
      
      return out_mat;
    }

} //end namespace fastmarching

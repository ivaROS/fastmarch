#include <fastmarch/dt_helper.h>

#include <opencv2/core/mat.hpp>

#include <fastmarch/signed_dist_test.h>
#include <opencv2/highgui.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>

//float INF = 1e20;

void print(std::vector<float> data, std::string name)
{
  std::cout << name << "\n[";
  
  for(int i = 0; i < data.size(); i++)
  {
    std::cout << data[i] << ",";
  }
  std::cout << "]\n\n";
}

std::vector<float> initialize(std::vector<float> data)
{
  std::vector<float> idata;
  for(int i = 0; i < data.size(); i++)
  {
    float v = (data[i]==0) ? INF : -INF;
    idata.push_back(v);
  }
  return idata;
}

void ssdt(std::vector<float>& f)
{
  int n = f.size();
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
}


//????
void adt(std::vector<float>& f)
{
  int n = f.size();
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
}



template <typename T>
class ArrayAccessor
{
public:
  virtual unsigned int size()=0;
  
	virtual T& operator[](unsigned int)=0;
};


template <typename T>
class MatRowArray : public ArrayAccessor<T>
{
  cv::Mat mat_;
  T* row_;
  
public:
  MatRowArray(cv::Mat mat, unsigned int row):
    mat_(mat)
  {
    row_ = (T*)mat_.ptr(row);
  }

  virtual unsigned int size()
  {
    return mat_.cols;
  }
  
	virtual T& operator[](unsigned int i)
  {
    return row_[i];
  }
};

template <typename T>
class MatColArray : public ArrayAccessor<T>
{
  cv::Mat mat_;
  unsigned int col_;
  
public:
  MatColArray(cv::Mat mat, unsigned int col):
    mat_(mat),
    col_(col)
  {
  }

  virtual unsigned int size()
  {
    return mat_.rows;
  }
  
	virtual T& operator[](unsigned int i)
  {
    return ((T*)mat_.ptr(i))[col_];
  }
};


template <typename T>
class MatAccessorInterface
{
  cv::Mat mat_;
  
public:
  MatAccessorInterface(cv::Mat mat):
    mat_(mat)
  {
  }
  
  MatRowArray<T> row(unsigned int i)
  {
    return MatRowArray<T>(mat_, i);
  }
  
  MatColArray<T> col(unsigned int i)
  {
    return MatColArray<T>(mat_, i);
  }
  
  cv::Mat getMat()
  {
    return mat_;
  }
};




//row-wise signed distance
void sdt(ArrayAccessor<float>& f)
{
  int n = f.size();
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
}



//row-wise signed distance
void dt_mine(ArrayAccessor<float>& f)
{
  int n = f.size();
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
}


/* dt of 1d function using squared distance */
void dt(ArrayAccessor<float>& f, ArrayAccessor<float>& d)
{
  int n = f.size();
  //float *d = new float[n];
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
}

void getSqrt(MatAccessorInterface<float>& in, MatAccessorInterface<float>& out)
{
  const cv::Mat& img = in.getMat();
//   if(img.isContinuous())
//   {
//     
//   }
//   else
  {

    for(int i = 0; i < img.rows; i++)
    {
      auto row_in = in.row(i);
      auto row_out = out.row(i);
      for(int j = 0; j < img.cols; j++)
      {
        row_out[j] = std::sqrt(row_in[j]);
      }
    }
  }
}


class DTApproach
{
public:
  virtual cv::Mat dt(cv::Mat image)=0;
  
  /* dt of 1d function using squared distance */
  static void dt(ArrayAccessor<float>& f, ArrayAccessor<float>& d)
  {
    int n = f.size();
    //float *d = new float[n];
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
  }
};

class InPlaceStandard : public DTApproach
{
public:
  cv::Mat dt(cv::Mat image) override
  {
    MatAccessorInterface<float> img(image);
    MatAccessorInterface<float> prepped(image.clone());

    //initialize distances of occupied points to 0 and the rest to a big number
    for(int i = 0; i < image.rows; i++)
    {
      auto row_in = img.row(i);
      auto row_out = prepped.row(i);
      for(int j = 0; j < image.cols; j++)
      {
        float v = (row_in[j]==0) ? INF : 0;
        row_out[j] = v;
      }
    }
    
    //Transform rows
    for(int i = 0; i < image.rows; i++)
    {
      auto row_in = prepped.row(i);
      DTApproach::dt(row_in, row_in);
    }
    
    //Transform columns
    for(int j = 0; j < image.cols; j++)
    {
      auto col_in = prepped.col(j);
      DTApproach::dt(col_in, col_in);
    }
    
    getSqrt(prepped,prepped);
    
    return prepped.getMat();
  }
};

class InPlaceReversed : public DTApproach
{
public:
  cv::Mat dt(cv::Mat image) override
  {
    MatAccessorInterface<float> img(image);
    MatAccessorInterface<float> prepped(image.clone());

    //initialize distances of occupied points to 0 and the rest to a big number
    for(int i = 0; i < image.rows; i++)
    {
      auto row_in = img.row(i);
      auto row_out = prepped.row(i);
      for(int j = 0; j < image.cols; j++)
      {
        float v = (row_in[j]==0) ? INF : 0;
        row_out[j] = v;
      }
    }
    
    for(int j = 0; j < image.cols; j++)
    {
      auto col_in = prepped.col(j);
      DTApproach::dt(col_in, col_in);
    }
    
    //initialize distances of occupied points to 0 and the rest to a big number
    for(int i = 0; i < image.rows; i++)
    {
      auto row_in = prepped.row(i);
      DTApproach::dt(row_in, row_in);
    }

    getSqrt(prepped,prepped);
    
    return prepped.getMat();
  }
};

class OutPlaceStandard : public DTApproach
{
public:
  cv::Mat dt(cv::Mat image) override
  {
    MatAccessorInterface<float> img(image);
    MatAccessorInterface<float> prepped(image.clone());

    //initialize distances of occupied points to 0 and the rest to a big number
    for(int i = 0; i < image.rows; i++)
    {
      auto row_in = img.row(i);
      auto row_out = prepped.row(i);
      for(int j = 0; j < image.cols; j++)
      {
        float v = (row_in[j]==0) ? INF : 0;
        row_out[j] = v;
      }
    }
    
    //Transform rows
    MatAccessorInterface<float> rows(image.clone());
    for(int i = 0; i < image.rows; i++)
    {
      auto row_in = prepped.row(i);
      auto row_out = rows.row(i);
      DTApproach::dt(row_in, row_out);
    }
    
    //Transform columns
    MatAccessorInterface<float> out(image.clone());
    for(int j = 0; j < image.cols; j++)
    {
      auto col_in = rows.col(j);
      auto col_out = out.col(j);
      DTApproach::dt(col_in, col_out);
    }
    
    MatAccessorInterface<float> sqrted(image.clone());
    getSqrt(out,sqrted);
    
    return sqrted.getMat();
  }
};

class DTHelperApproach : public DTApproach
{
  fastmarching::DTHelper helper_;
  
public:
  cv::Mat dt(cv::Mat image) override
  {
    return helper_.run(image);
  }

  
};

cv::Mat accessorOPDT(cv::Mat image)
{
  MatAccessorInterface<float> img(image);
  MatAccessorInterface<float> prepped(image.clone());

  //initialize distances of occupied points to 0 and the rest to a big number
  for(int i = 0; i < image.rows; i++)
  {
    auto row_in = img.row(i);
    auto row_out = prepped.row(i);
    for(int j = 0; j < image.cols; j++)
    {
      float v = (row_in[j]==0) ? INF : 0;
      row_out[j] = v;
    }
  }
  
  MatAccessorInterface<float> rows(image.clone());
  //initialize distances of occupied points to 0 and the rest to a big number
  for(int i = 0; i < image.rows; i++)
  {
    auto row_in = prepped.row(i);
    auto row_out = rows.row(i);
    dt(row_in, row_out);
  }
  
  MatAccessorInterface<float> out(image.clone());
  for(int j = 0; j < image.cols; j++)
  {
    auto col_in = rows.col(j);
    auto col_out = out.col(j);
    dt(col_in, col_out);
  }
  
  return out.getMat();
}


void test0()
{
  std::vector<float> data;
  int n = 20;
  
  for(int i = 0; i < n; i++)
  {
    data.push_back(0);
  }
  data[0] = data[6] = data[7] = data[8] = data[9] = data[10] = data[15] = 1;
  
  print(data, "original");
  
  auto idata = initialize(data);
  print(idata, "initialized");
  
  sdt(idata.data(), idata.size());
  //ssdt(idata);

  print(idata, "sdt");
 
}


class ApproachTester
{
  cv::Mat result_;
  std::string name_;
  std::chrono::duration< double, std::milli > ttime_;
  
protected:
  DTApproach* approach_;

public:

  ApproachTester(std::string name):
    name_(name)
  {}
  
  void runTest(cv::Mat image)
  {
    auto t1 = std::chrono::high_resolution_clock::now();
    result_ = approach_->dt(image);
    auto t2 = std::chrono::high_resolution_clock::now();
    ttime_ = t2 - t1;
  }
  
  void printTime()
  {
    std::cout << name_ << ": " << ttime_.count() << "ms\n";
  }
  
  void showResult()
  {
    cv::Mat result_viz;
    cv::normalize(result_, result_viz, 1, 0, cv::NORM_INF);
    cv::imshow(name_ + ":result", result_viz);
  }
};

template <typename T>
class TypedApproachTester : public ApproachTester
{
  T approach_impl_;
  
public:
  TypedApproachTester(std::string name):
    ApproachTester(name)
  {
    approach_ = &approach_impl_;
  }
};

void compareInOutPlace(cv::Mat image)
{
  TypedApproachTester<InPlaceStandard> ipt("in-place");
  TypedApproachTester<OutPlaceStandard> opt("out-of-place");
  TypedApproachTester<DTHelperApproach> dth("first-try");
  TypedApproachTester<InPlaceReversed> iprt("in-place-reversed");

  std::vector<ApproachTester*> tests = {&ipt, &opt, &dth, &iprt};
  
  for(int i = 0; i < 20; i++)
  {
    for(auto test : tests)
    {
      test->runTest(image);
    }
    
    for(auto test : tests)
    {
      test->printTime();
    }
  }
  
  for(auto test : tests)
  {
    test->showResult();
  }

  //InPlaceStandard ip;
  //OutPlaceStandard op;
  
  //std::vector<DTApproach*> approaches = {&ip, &op};
  
  //auto t1 = std::chrono::high_resolution_clock::now();
  //ApproachTester ipt("in-place");
  //ipt.runTest(ip, image);

  //auto t2 = std::chrono::high_resolution_clock::now();
  //ApproachTester opt("out-of-place");
  //opt.runTest(op, image);
  
  //auto t3 = std::chrono::high_resolution_clock::now();

  //Calculate elapsed time for the computations
  //std::chrono::duration<double, std::milli> iptime = t2 - t1;
  //std::chrono::duration<double, std::milli> optime = t3 - t2;

  //std::cout << "In-place computation: " << iptime.count() << "ms; Out-of-place computation: " << optime.count() << "ms\n";
  //ipt.showResult();
  //opt.showResult();
}



void test1()
{
  cv::Mat input = cv::Mat::zeros(128,512,CV_32FC1);
  input.at<float>(64,256) = 1;
  compareInOutPlace(input);
//   cv::Mat result = accessorOPDT(input);
//   cv::Mat result_viz;
//   cv::normalize(result, result_viz, 1, 0, cv::NORM_INF);
//   
   cv::imshow("input", input);
//   cv::imshow("result", result_viz);
}



int main(void)
{
  test0();
  test1();
  
  cv::waitKey();

  return 0;
}

#ifndef CAMEAR_H
#define CAMERA_H
#include <string>
#include <opencv2/opencv.hpp>
using namespace std;

struct Frame {
  int     frame_id;
  cv::Mat data;
};	

class Camera {
 public:
  Camera();
  ~Camera();

  void SetCameraParameters(double dx, double dy, int nRows, int nCols, string model);
  
  double dx_;
  double dy_;
  int nRows_;
  int nCols_;
  string model_;
  bool dx_init;
  bool dy_init;
  bool nRows_init;
  bool nCols_init;
    
};


#endif
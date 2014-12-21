#ifndef MONOSLAM_H
#define MONOSLAM_H

#include <cmath>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <string>
#include <Eigen/Eigen>
#include <iostream>
#include "vars.h"
#include "FilterBank.h"
#include "KalmanFilter.h"
#include "Camera.h"
#include "MotionModel.h"


using namespace std;
using namespace Eigen;

class MapManagement;
class Kalman;
class DataAssociator;

class MonoSLAM{

public:
	MonoSLAM();
	~MonoSLAM();

	void Init(const string &config_path, double* dx, double* dy,
		int* nRows, int* nCols, string* model, int init_frame_id, int init_max);
	void Init_filterbank(FilterBank *KalmanFilterBank_);
	void GoOneStep(int step, int init_frame_id, int init_max, Frame frame_last, Frame frame_next);

	struct feature_info
	{
		cv::Mat patch_when_initialized;
		MatrixXd patch_when_matching;
		VectorXd r_wc_when_initialized;
		MatrixXd R_wc_when_initialized;
		MatrixXd uv_when_initialized;
		int half_patch_size_when_initialized;
		int half_patch_size_when_matching;
		int times_predicted;
		int times_measured;
		int init_frame;
    int ind_feat;
		VectorXd init_measurement;
		string type;
		VectorXd yi;
	};
	vector<feature_info> features_info;
	
	// Feature points measured for every image frame
	struct uv_initialized
	{
		int index_frame;
		MatrixXd initialized_features_per_frame;
		int features_size;
		bool isempty;
	};
	vector<uv_initialized> initialized_features;

	bool input_mode;
	int min_number_of_features_in_image;
	double A, B;
	Matrix <double, Dynamic, 2> measurements;
	Matrix <double, Dynamic, 2> predicted_measurements;
  // To compare the results between scopis calibration and this calibration
  Matrix <double, Dynamic, 2> error_coordinates; // in pixel
  MatrixXd S_predicted_;
	VectorXd mu, likelihood_ratio;
	// Final state and convariance output value
	VectorXd x_k_k_output;
	MatrixXd p_k_k_output;
	string input_name;

	MapManagement            *MapManager;
	FilterBank               *KalmanFilterBank_;
	Camera                   *Cam;
	Kalman                   *Kalman_;
	MotionModel              *MotionModel_;
	DataAssociator           *DataAssociator_;

};

#endif
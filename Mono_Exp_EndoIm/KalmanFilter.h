#ifndef KALMANFILTER_H
#define KALMANFILTER_H
#include<string>
#include<Eigen/Eigen>

using namespace std;
using namespace Eigen;

class KalmanFilter 
{
public: 
	KalmanFilter();
	~KalmanFilter();
	void InitFilterFields(double std_a,
	                     double std_alpha,
	                     double standard_deviation_z,
						 VectorXd x_k_km1,
						 MatrixXd p_k_km1,
						 string type);
	string type;

	// State vector and covariance matrix
	VectorXd x_k_k_;
	MatrixXd p_k_k_;
	double std_a_;
	double std_alpha_;
	double std_z_;
	//double standard_deviation_z_;
	// Number of calibration, pisition and orientation(cartesian and quaternion) parameters
	int num_calibration, num_t_wc, num_q_wc, num_r_wc, num_cam_para;
	// Files below not sure about matrix or vector
	VectorXd x_k_km1_;
	MatrixXd p_k_km1_;
	Matrix <double, Dynamic, 2> predicted_measurements;
	size_t num_measurements;
	MatrixXd H_predicted_;
	MatrixXd R_predicted_;
	MatrixXd S_predicted_;
	size_t Rownum_S_predicted;
	size_t Colnum_S_predicted;
	MatrixXd S_matching_;
	MatrixXd z_;
	MatrixXd h_;
	MatrixXd H_matching_;
	MatrixXd measurements_;
	MatrixXd R_matching_;
	MatrixXd x_k_k_mixing_estimate_;
	MatrixXd p_k_k_mixing_covariance_;

};

#endif
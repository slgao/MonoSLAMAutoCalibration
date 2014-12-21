#ifndef MAPMANAGEMENT_H
#define MAPMANAGEMENT_H

#include "Monoslam.h"
#include "MathUtil.h"

using namespace std;
using namespace Eigen;


class MapManagement{
public: 
	MapManagement();
	~MapManagement();


	void map_management_filter_bank(int step, Frame frame, MonoSLAM *mono_slam);
	void update_features_info(MonoSLAM* mono_slam);
    void inversedepth_2_cartesian_filter_bank(MonoSLAM* mono_slam);
	Vector3d inversedepth2cartesian(VectorXd inverse_depth);
    void delete_features_filter_bank(MonoSLAM* mono_slam);
	
	void delete_a_feature(VectorXd& x_k_k, MatrixXd& p_k_k, int FeatToDelete, vector<MonoSLAM::feature_info> features_info);


    void initialize_features_filter_bank(int step, int num_feat_to_init, MonoSLAM *mono_slam, Frame frame);
	void initialize_a_feature_other_filter_bank(int step, Frame frame, MonoSLAM *mono_slam, cv::KeyPoint *feature_TBadded);
	void add_features_inverse_depth(cv::KeyPoint feature, VectorXd *x_k_k, MatrixXd *p_k_k, 
		Camera *cam, double std_z, double initial_rho, double std_rho, VectorXd *newFeature);
	void hinv(cv::KeyPoint feature, VectorXd Xv, Camera * cam, double initial_rho, VectorXd *newFeature);
	// Undistort the detected points
	void undistort_fm(cv::KeyPoint *undistort_feature, cv::KeyPoint feature, Camera *cam, VectorXd calibration);
	MatrixXd add_a_feature_covariance_inverse_depth(MatrixXd P, cv::KeyPoint feature, VectorXd Xv, double std_z, double std_rho, Camera *cam);
	MatrixXd jacob_undistort_fm(Camera *cam, cv::KeyPoint feature, VectorXd calibration);

	
	void find_matched_measurements(int *num_of_matched, MatrixXd measurements);
    MatrixXd predict_camera_measurements(VectorXd x_k_k, KalmanFilter *kalmanfilter, MonoSLAM *mono_slam);
	// Convert from quaternion to rotation matrix in cartesian representation
	bool hi_cartesian(VectorXd *zi, VectorXd yi3d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
						vector<MonoSLAM::feature_info> features_info, VectorXd cali_para);
	bool hi_inverse_depth(VectorXd *zi, VectorXd yi6d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
						vector<MonoSLAM::feature_info> features_info, VectorXd cali_para);
	VectorXd hu(VectorXd yi, Camera *cam, VectorXd cali_para);
	// Project the undistorted image coordinates to the distorted ones
	Vector2d distort_fm(VectorXd uv, Camera *cam, VectorXd cali_para);

	const int max_attempts, half_patch_size_when_initialized, 
		half_patch_size_when_matching, max_initialization_attempts, 
		horizontal_div, vertical_div;
	int horizontal_size, vertical_size;
	bool imsize_determined, feature_added;
	MatrixXd predicted_measurements;
	
};

#endif
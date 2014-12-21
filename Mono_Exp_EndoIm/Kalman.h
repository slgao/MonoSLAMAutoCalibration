#ifndef KALMAN_H
#define KALMAN_H
#include "Monoslam.h"
#include "MathUtil.h"

class Kalman
{
public:
	Kalman();
	~Kalman();

	// The function predict_camera_measurements and subfunctions are also defined in MapManagementFilterBank file.
	// Structure of the software needs to be reorgnized afterwards for optimization. One idea is to make a predict camera measurements class and 
	// this class should not include Monoslam.h file and should not have MonoSLAM type input argument, use the class in Kalman and MapManagementFilterBank file 
	// in the end.
    // ------------------------------------------------------------------------------------------------------------
    void predict_camera_measurements(VectorXd x_k_k, KalmanFilter *kalmanfilter, MonoSLAM *mono_slam, Matrix <double, Dynamic, 2> *predicted_measurements);
	bool hi_cartesian(VectorXd *zi, VectorXd yi3d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
						vector<MonoSLAM::feature_info> features_info, VectorXd cali_para);
	bool hi_inverse_depth(VectorXd *zi, VectorXd yi6d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
						vector<MonoSLAM::feature_info> features_info, VectorXd cali_para);
  bool hi_cartesian_scopis(VectorXd *zi, VectorXd yi3d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
						vector<MonoSLAM::feature_info> features_info, VectorXd cali_para);
	bool hi_inverse_depth_scopis(VectorXd *zi, VectorXd yi6d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
						vector<MonoSLAM::feature_info> features_info, VectorXd cali_para);
	VectorXd hu(VectorXd yi, Camera *cam, VectorXd cali_para);
	// Project the undistorted image coordinates to the distorted ones
	Vector2d distort_fm(VectorXd uv, Camera *cam, VectorXd cali_para);
    // ------------------------------------------------------------------------------------------------------------


	void CalculateDerivatives(VectorXd x_k_km1, MatrixXd PredictedMeasurements, Camera *cam, vector<MonoSLAM::feature_info> features_info, MatrixXd *H_predicted);
    
	// Calculation of Hi derivative in inverse depth
	void CalculateHiInverseDepth(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, size_t FeatIndex, 
		vector<MonoSLAM::feature_info> features_info, MatrixXd *H_predicted);
	void dhd_df(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex);
	MatrixXd dhd_dhu(Camera *cam, VectorXd zi_d, VectorXd calibration);
	void JacobUndistortFM(Camera *cam, VectorXd uvd, VectorXd calibration, MatrixXd *J_undistort);
	VectorXd dhu_df(Camera *cam, VectorXd Xv_km1_k, VectorXd yi);
	void dhd_dCx(MatrixXd *H_predicted, size_t FeatIndex);
	void dhd_dCy(MatrixXd *H_predicted, size_t FeatIndex);

	// This function is not consistent with the theory explanationin the Javier book P142. A.101
	void dhd_dk1(Camera *cam, cv::KeyPoint feature, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex);
	// This function is not consistent with the theory explanationin the Javier book P142. A.101
	void dhd_dk2(Camera *cam, cv::KeyPoint feature, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex);

	void undistort_fm(cv::KeyPoint *undistort_feature, cv::KeyPoint feature, Camera *cam, VectorXd calibration);
	void dh_dxv(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t FeatIndex);
	MatrixXd dh_drw(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi);
	MatrixXd dh_dhrl(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi);
	MatrixXd dhu_dhrl(Camera *cam, VectorXd Xv_km1_k, VectorXd yi);
	MatrixXd dhrl_drw(VectorXd Xv_km1_k, VectorXd yi);
	MatrixXd dh_dqwr(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi);
	MatrixXd dhrl_dqwr(VectorXd Xv_km1_k, VectorXd yi);
	VectorXd QConj(VectorXd q);

	void dh_dy(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t index_insertion, size_t FeatIndex);
	MatrixXd dhrl_dy(VectorXd Xv_km1_k, VectorXd yi);
	// Calculation of Hi derivative in inverse depth

	// Calculation of Hi derivative in Cartesian
	void CalculateHiCartesian(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, size_t FeatIndex, 
	vector<MonoSLAM::feature_info> features_info, MatrixXd *H_predicted);

	// Calculation of dhd_df_xyz
	void dhd_df_xyz(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex);
	VectorXd dhu_df_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi);
	// Calculation of dhd_df_xyz

	// Calculation of dh_drw_xyz
	void dh_dxv_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t FeatIndex);
	MatrixXd dh_drw_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi);
	MatrixXd dh_dhrl_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi);
	MatrixXd dhu_dhrl_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi);
	MatrixXd dhrl_drw_xyz(VectorXd Xv_km1_k);
	MatrixXd dh_dqwr_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi);
	MatrixXd dhrl_dqwr_xyz(VectorXd Xv_km1_k, VectorXd yi);
	// Calculation of dh_drw_xyz

	// Calculation of dh_dy_xyz
	void dh_dy_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t index_insertion, size_t FeatIndex);
	MatrixXd dhrl_dy_xyz(VectorXd Xv_km1_k);
	// Calculation of dh_dy_xyz
	// Calculation of Hi derivative in Cartesian

	// Extended kalman filter prediction step
	void EKF_Prediction(MonoSLAM *mono_slam);
	void EKF_Update(vector<KalmanFilter *> FilterBank);
	void Update(VectorXd x_km1_k, MatrixXd p_km1_k, MatrixXd H, MatrixXd R, VectorXd z, VectorXd h, VectorXd *x_k_k, MatrixXd *p_k_k);
	MatrixXd NormJac(VectorXd q);


};

#endif
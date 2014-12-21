#include "Kalmanfilter.h"

KalmanFilter::KalmanFilter():num_calibration(5),
						     num_t_wc(3), num_q_wc(4), num_r_wc(3), num_cam_para(18)// number of prameters initialized to 5, 3, 4, 3
{

}

KalmanFilter::~KalmanFilter()
{



}

void KalmanFilter::InitFilterFields(double std_a,
	                                double std_alpha,
	                                double standard_deviation_z,
									VectorXd x_k_k,
									MatrixXd p_k_k,
									string type)
{
	std_a_ = std_a;
	std_alpha_ = std_alpha;
	std_z_ = standard_deviation_z;
	x_k_k_ = x_k_k;
	p_k_k_ = p_k_k;
	this->type = type;
}

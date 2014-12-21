#include "FilterBank.h"
#include <math.h>

#define eps 2.220446049250313e-016

FilterBank::FilterBank():f_start(100), f_end(610), k1_start(0.02), k1_end(0.06),
	k2_start(0.003), k2_end(0.015), standard_deviation_z(1.0), std_a(0.005), std_alpha(0.005),
	f_deviation(30), k1_deviation(0.04), k2_deviation(0.006), Cx_0(160), Cy_0(120), type("constant_velocity")
  // Defaut
  //f_start(100), f_end(610), k1_start(0.02), k1_end(0.06),
	//k2_start(0.003), k2_end(0.015), standard_deviation_z(1.0), std_a(0.005), std_alpha(0.005),
	//f_deviation(30), k1_deviation(0.04), k2_deviation(0.006), Cx_0(160), Cy_0(120), type("constant_velocity")
{
	x_k_k.resize(18);
	p_k_k.resize(18,18);
	filter = NULL;
}

FilterBank::~FilterBank()
{
	if (filter != NULL)
		delete  filter;
	while (!FilterBank_.empty())
		FilterBank_.pop_back();
}

void FilterBank::initialize_x_and_p()
{
	// Initial velocity and calibration parameters values
  // Defaut:
  //double v_0 = 0; double std_v_0 = 0; double w_0 = 1e-15; double std_w_0 = 0;
	//double f_0 = 210; double std_f_0 = 7.5; double Cx_0 = 160; double std_Cx_0 = 3.3;
	//double Cy_0 = 120; double std_Cy_0 = 3.3; double k1_0 = 0.069; double std_k1_0 = 0.01;
	//double k2_0 = 0.011; double std_k2_0 = 0.0015;

	double v_0 = eps; double std_v_0 = 0; double w_0 = 1e-15; double std_w_0 = 0;
  //                                            1e-15
	double f_0 = 210; double std_f_0 = 7.5; double Cx_0 = 360; double std_Cx_0 = 3.3;
	double Cy_0 = 290; double std_Cy_0 = 3.3; double k1_0 = 0.069; double std_k1_0 = 0.01;
	double k2_0 = 0.011; double std_k2_0 = 0.0015;

	// Initialize state vector and covariance matrix
	x_k_k << f_0, Cx_0, Cy_0, k1_0, k2_0, 0, 0, 0, 1, 0, 0, 0, v_0, v_0, v_0, w_0, w_0, w_0;
	p_k_k.setZero();
	p_k_k(0,0)=pow(std_f_0, 2);
	p_k_k(1,1)=pow(std_Cx_0, 2);
	p_k_k(2,2)=pow(std_Cy_0, 2);   
	p_k_k(3,3)=pow(std_k1_0, 2);   
	p_k_k(4,4)=pow(std_k2_0, 2);   
	p_k_k(5,5)=eps;
	p_k_k(6,6)=eps;
	p_k_k(7,7)=eps;
	p_k_k(8,8)=eps;
	p_k_k(9,9)=eps;
	p_k_k(10,10)=eps;
	p_k_k(11,11)=eps;
	p_k_k(12,12)=pow(std_v_0, 2);
	p_k_k(13,13)=pow(std_v_0, 2);
	p_k_k(14,14)=pow(std_v_0, 2);
	p_k_k(15,15)=pow(std_w_0, 2);
	p_k_k(16,16)=pow(std_w_0, 2);
	p_k_k(17,17)=pow(std_w_0, 2);
}

void FilterBank::initialize_filterbank()
{
	int f_ceil = (int)ceil((f_end - f_start)/f_deviation);
	int k1_ceil = (int)ceil((k1_end - k1_start)/k1_deviation);
	int k2_ceil = (int)ceil((k2_end - k2_start)/k2_deviation);

	//values combination of focal length and two distortion parameters to initialize filter bank
	for (int i_f = 0; i_f <= f_ceil; i_f ++){
		for (int i_k1 = 0; i_k1 <= k1_ceil; i_k1 ++){
			for (int i_k2 = 0; i_k2 <= k2_ceil; i_k2 ++){
				x_k_k(0) = f_start + i_f * f_deviation;
				x_k_k(3) = k1_start + i_k1 * k1_deviation; x_k_k(4) = k2_start + i_k2 * k2_deviation;
				filter = new KalmanFilter();
				filter->InitFilterFields(std_a, std_alpha, standard_deviation_z, x_k_k, p_k_k, type);
				FilterBank_.push_back(filter);
			}	  

		}
	}
	filter_size = FilterBank_.size();

}




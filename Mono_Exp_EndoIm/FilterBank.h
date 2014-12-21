#ifndef FILTERBANK_H
#define FILTERBANK_H
#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include "KalmanFilter.h"

using namespace Eigen;
using namespace std;

class FilterBank {
public:
	FilterBank();
	~FilterBank();
	void initialize_x_and_p();
	void initialize_filterbank();
	
	// state vector and covariance vector for initialization
	VectorXd x_k_k; 
	MatrixXd p_k_k;
	double standard_deviation_z;
	double f_start;
	double f_end;
	double f_deviation;
	double k1_start;
	double k1_end;
	double k1_deviation;
	double k2_start;
	double k2_end;
	double k2_deviation;
	double std_a;
	double std_alpha;
    double Cx_0;
	double Cy_0;
	string type;
	size_t filter_size;

	KalmanFilter *filter;
	vector<KalmanFilter*> FilterBank_;

};

#endif
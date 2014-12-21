#ifndef MOTIONMODEL_H
#define MOTIONMODEL_H
#include <math.h>
#include <vector>
#include <string>
#include <Eigen/Eigen>

#define eps 2.220446049250313e-016
using namespace Eigen;
using namespace std;


class MotionModel {
public:

	MotionModel();
	~MotionModel();
	void predict_state_and_covariance(VectorXd x_k_k, MatrixXd p_k_k, string type, double std_a,
		double std_alpha, VectorXd *X_km1_k, MatrixXd *P_km1_k);
	void func_Q(VectorXd Xv, MatrixXd Pn, double delta_t, string type, MatrixXd *G, MatrixXd *Q);

	// The three functions below are probably not used
	// Convert unit-quaternion to homogeneous transform
	MatrixXd q2tr(VectorXd q);
	// Convert a homogeneous transform matrix to roll/pitch/yaw angles
	VectorXd tr2rpy(MatrixXd m);
	void dq_by_deuler(VectorXd euler_angles, MatrixXd *G);
	// ---------------------------------------------------------------

	void fv(VectorXd x_k_k, double delta_t, string type, double std_a, double std_alpha, VectorXd *Xv_km1_k);
	VectorXd V2Q(VectorXd V);
	VectorXd Qprod(VectorXd Q, VectorXd P);
	void dfv_by_dxv(VectorXd Xv, double delta_t, string type, MatrixXd *dfv_by_dxvRES);
	MatrixXd dq3_by_dq2(VectorXd q1_in);
	MatrixXd dq3_by_dq1(VectorXd q2_in);
	MatrixXd dqomegadt_by_domega(VectorXd omega, double delta_t);
	double dq0_by_domegaA(double omegaA, double omega_norm, double delta_t);
	double dqA_by_domegaA(double omegaA, double omega_norm, double delta_t);
	double dqA_by_domegaB(double omegaA, double omegaB, double omega_norm, double delta_t);

	
};


#endif 

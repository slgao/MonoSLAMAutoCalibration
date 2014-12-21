#include "MotionModel.h"



MotionModel::MotionModel()
{


}

MotionModel::~MotionModel()
{


}

void MotionModel::predict_state_and_covariance(VectorXd x_k_k, MatrixXd p_k_k, string type, double std_a,
	double std_alpha, VectorXd *X_km1_k, MatrixXd *P_km1_k)
{
	size_t x_k_size = x_k_k.size(), p_k_size = p_k_k.rows();
	double delta_t = 1, linear_acceleration_noise_covariance, angular_acceleration_noise_covariance;
	VectorXd Xv_km1_k(18), Pn_diag(6);
	MatrixXd F = MatrixXd::Identity(18,18), // The camera state size is assumed to be 18
		Pn, Q, G;

	// Camera motion prediction
	fv(x_k_k.head(18), delta_t, type, std_a, std_alpha, &Xv_km1_k);

	// Feature prediction
	(*X_km1_k).resize(x_k_size);
	// Add the feature points to the state vector after cam motion prediction
	(*X_km1_k) << Xv_km1_k, x_k_k.tail(x_k_size - 18);

	// State transition equation derivatives
	dfv_by_dxv(x_k_k.head(18), delta_t, type, &F);

	// State noise
	linear_acceleration_noise_covariance = (std_a * delta_t) * (std_a * delta_t);
	angular_acceleration_noise_covariance = (std_alpha * delta_t) * (std_alpha * delta_t);
	Pn_diag << linear_acceleration_noise_covariance, linear_acceleration_noise_covariance, linear_acceleration_noise_covariance, 
		angular_acceleration_noise_covariance, angular_acceleration_noise_covariance, angular_acceleration_noise_covariance;
	Pn = Pn_diag.asDiagonal();
	func_Q(x_k_k.head(18), Pn, delta_t, type, &G, &Q);

	// P_km1_k resize and update 
	// With the property of symmetry for the covariance matrix, only the Pxx and Pxy 
	// which means the camera state covariance and camera feature points covariance can be calculated
	(*P_km1_k).resize(p_k_size,p_k_size);
	(*P_km1_k).block(0,0,18,18) = F * p_k_k.block(0,0,18,18) * F.transpose() + Q;
	(*P_km1_k).block(0,18,18,p_k_size - 18) = F * p_k_k.block(0,18,18,p_k_size - 18);
	(*P_km1_k).block(18,0,p_k_size - 18,18) = p_k_k.block(18,0,p_k_size - 18,18) * F.transpose();
	(*P_km1_k).block(18,18,p_k_size - 18,p_k_size - 18) = p_k_k.block(18,18,p_k_size - 18,p_k_size - 18);

}

void MotionModel::func_Q(VectorXd x_k_k, MatrixXd Pn, double delta_t, string type, MatrixXd *G, MatrixXd *Q)
{
	// This function calculates the noise Q, see javier book P128. P129. A.11
	// input: x_k_k.head(18), Pn, delata_t, type
	// output: G, Q
	// Pn is the input noise covariance, G is the derivatives of the dynamic model by the 
	// input parameters.
	double std_f_0 = 0, std_Cx_0 = 0, std_Cy_0 = 0, std_k1_0 = 0, std_k2_0 = 0;
	VectorXd qOld, omegaOld;
	if (type.compare("constant_position_and_orientation_location_noise") == 0)
	{
		qOld = x_k_k.segment(3,4);
		(*G) = MatrixXd::Zero(13,6);
		(*G).block(0,0,3,3) = MatrixXd::Identity(3,3) * delta_t;
		dq_by_deuler(tr2rpy(q2tr(qOld)), G);
	}
	else
	{
		omegaOld = x_k_k.segment(15,3);
		qOld = x_k_k.segment(8,4);
		(*G) = MatrixXd::Zero(13,6);
		(*G).block(7,0,3,3) = MatrixXd::Identity(3,3);
		(*G).block(10,3,3,3) = MatrixXd::Identity(3,3);
		(*G).block(0,0,3,3) = MatrixXd::Identity(3,3) * delta_t;
		(*G).block(3,3,4,3) = dq3_by_dq1(qOld) * dqomegadt_by_domega(omegaOld, delta_t); // P129. A.11
	}

	(*Q) = MatrixXd::Zero(18,18);
	(*Q)(0,0) = std_f_0 * std_f_0;
	(*Q)(1,1) = std_Cx_0 * std_Cx_0;
	(*Q)(2,2) = std_Cy_0 * std_Cy_0;
	(*Q)(3,3) = std_k1_0 * std_k1_0;
	(*Q)(4,4) = std_k2_0 * std_k2_0;
	(*Q).block(5,5,13,13) = (*G) * Pn * (*G).transpose(); // P128. A.4

}

// The three functions below are probaby not used
MatrixXd MotionModel::q2tr(VectorXd q)
{
	double s, x, y, z;
	Matrix3d r;
	MatrixXd t(4,4);
	s = q(0);
	x = q(1);
	y = q(2);
	z = q(3);
	r << 1 - 2 * (y * y + z * z), 2 * (x * y - s * z), 2 * (x * z + s * y),
		2 * (x * y + s * z), 1 - 2 * (x * x + z * z), 2 * (y * z - s * x),
		2 * (x * z - s * y), 2 * (y * z + s * x), 1 - 2 * (x * x + y * y);
	t = MatrixXd::Identity(4,4);
	t.block(0,0,3,3) = r;
	return t;
}

VectorXd MotionModel::tr2rpy(MatrixXd m)
{
	VectorXd rpy(3);
	if (abs(m(0,0)) < eps && abs(m(1,0)) < eps)
		rpy << 0, atan2(-m(2,0), m(0,0)), atan2(-m(1,2), m(1,1));
	else
	{
		double sp, cp;
		sp = sin(atan2(m(1,0), m(0,0)));
		cp = cos(atan2(m(1,0), m(0,0)));
		rpy << atan2(m(1,0), m(0,0)), atan2(-m(2,0), cp * m(0,0) + sp * m(1,0)), 
			atan2(sp * m(0,2) - cp * m(1,2), cp * m(1,1) - sp * m(0,1));
	}
	return rpy;
}

void MotionModel::dq_by_deuler(VectorXd euler_angles, MatrixXd *G)
{
	double phi, theta, psi;
	phi = euler_angles(0);
	theta = euler_angles(1);
	psi = euler_angles(2);

	(*G).block(3,3,4,3) << (0.5)*(-sin(phi/2)+cos(phi/2)), (0.5)*(-sin(theta/2)+cos(theta/2)), (0.5)*(-sin(psi/2)+cos(psi/2)),
		(0.5)*(+cos(phi/2)+sin(phi/2)), (0.5)*(-sin(theta/2)-cos(theta/2)), (0.5)*(-sin(psi/2)-cos(psi/2)), 
		(0.5)*(-sin(phi/2)+cos(phi/2)), (0.5)*(+cos(theta/2)-sin(theta/2)), (0.5)*(-sin(psi/2)+cos(psi/2)), 
		(0.5)*(-sin(phi/2)-cos(phi/2)), (0.5)*(-sin(theta/2)-cos(theta/2)), (0.5)*(+cos(psi/2)+sin(psi/2));
}

void MotionModel::fv(VectorXd x_k_k, double delta_t, string type, double std_a, double std_alpha, VectorXd *Xv_km1_k)
{
	// This function updates the state for the camera. See Javier Civera book P39. 3.2.1 or P129 A.9
	// It assumes the camera has a constant linear and angular velocity, so the 
	// linear and angular accelerations are zero.

	VectorXd calibration, rW, qWR, vW, wW, Q, QP; // The camera state size is assumed to be 18
	calibration = x_k_k.head(5); // Five intrinsic parameters are assumed
	rW = x_k_k.segment(5,3); // Camera optical center position
	qWR = x_k_k.segment(8,4); // Quaternion defines orientation
	vW = x_k_k.segment(12,3); // Linear velocity
	wW = x_k_k.segment(15,3); // Angular velocity

	if (type.compare("constant_velocity") == 0)
	{
		Q = V2Q(wW * delta_t); // Converts rotation vector to quaternion
		QP = Qprod(qWR, Q); // Cross product of qWR and q((wW+OmegaW)*delta_t), where OmegaW is set to 0
		*Xv_km1_k << calibration, rW + vW * delta_t, QP, vW, wW; // State update for the camera
	}

	// The types below are probably not used
	else if (type.compare("constant_orientation") == 0)
	{

	}
	else if (type.compare("constant_position") == 0)
	{

	}
	else if (type.compare("constant_position_and_orientation") == 0)
	{

	}
	else if (type.compare("constant_position_and_orientation_location_noise") == 0)
	{

	}
}

VectorXd MotionModel::V2Q(VectorXd V)
{
	// Converts from rotation vector to quaternion representation
	// Javier Civera book P130. A.15
	size_t V_size;
	double norm = 0, v_norm = 0;
	VectorXd Q(4), V_norm;
	V_size = V.size();
	// Norm of the V
	//for (i = 0; i < V_size; ++i)
	//{
	//	norm += V(i) * V(i);
	//}
	//norm = sqrt(norm);

	norm = V.norm();
	if (norm < eps)
		Q << 1, 0, 0, 0;
	else
	{
		V_norm = V / norm;
		//for (i = 0; i < V_size; ++i)
		//{
		//	v_norm += V_norm(i) * V_norm(i);
		//}
		//v_norm = sqrt(v_norm);

		v_norm = V_norm.norm();
		// Quaternion cos(theta/2)+sin(theta/2)uxI+sin(theta/2)uyJ+sin(theta/2)uzK
		Q << cos(norm / 2), sin(norm / 2) * V_norm / v_norm;
	}
	return Q; // Q is a 4x1 vector
}

VectorXd MotionModel::Qprod(VectorXd Q, VectorXd P)
{
	double a, x;
	VectorXd QP(4), QP2;
	Vector3d u, v;
	a = Q(0); x = P(0);
	v = Q.segment(1,3);
	u = P.segment(1,3);
	QP2 = a * u + x * v + v.cross(u); // cross of v and u
	QP << a * x - v.transpose() * u, QP2;
	return QP; // QP should be a 4x1 vector
}

void MotionModel::dfv_by_dxv(VectorXd x_k_k, double delta_t, string type, MatrixXd *dfv_by_dxvRES)
{
	// This function calculates the F_k which is the Jacobian of f_k wrt. the state vector x_k_k-1
	// Javier Civera book P72. 4.3.1 or P129. A.10
	// Input : x_k_k.head(18), delta_t, type
	// Output: F_k

	VectorXd omegaOld, qOld, qwt;
	omegaOld = x_k_k.segment(15,3); // Angular velocity
	qOld = x_k_k.segment(8,4); // Orientation quaternion representation

	qwt = V2Q(omegaOld * delta_t); // Convert omegaOld to quaternion representation
	(*dfv_by_dxvRES).block(8,8,4,4) = dq3_by_dq2(qwt); // dq_k+1^WC/dq_k^WC

	if (type.compare("constant_velocity") == 0)
	{
		(*dfv_by_dxvRES).block(5,12,3,3) = MatrixXd::Identity(3,3) * delta_t; // Delta_t * Identity
		(*dfv_by_dxvRES).block(8,15,4,3) = dq3_by_dq1(qOld) * dqomegadt_by_domega(omegaOld, delta_t); // dq_k+1^WC/dw_k^C chain rule
	}
	/*
	else if (type.compare("constant_orientation") == 0)
	{
	dfv_by_dxvRES.block(8,15,4,3) = MatrixXd::Zero(4,3);
	dfv_by_dxvRES.block(15,15,3,3) = MatrixXd::Zero(3,3);
	}

	else if (type.compare("constant_position") == 0)
	{
	dfv_by_dxvRES.block(5,12,3,3) = MatrixXd::Zero(3,3);
	dfv_by_dxvRES.block(12,12,3,3) = MatrixXd::Zero(3,3);
	}

	else if (type.compare("constant_position_and_orientation") == 0)
	{
	dfv_by_dxvRES.block(8,15,4,3) = MatrixXd::Zero(4,3);
	dfv_by_dxvRES.block(5,12,3,3) = MatrixXd::Zero(3,3);
	dfv_by_dxvRES.block(15,15,3,3) = MatrixXd::Zero(3,3);
	dfv_by_dxvRES.block(12,12,3,3) = MatrixXd::Zero(3,3);
	}
	*/
}

MatrixXd MotionModel::dq3_by_dq2(VectorXd q1_in)
{
	// This function calculates q_k+1^WC/q_k^WC
	// q1_in = q((w_k^C + omega^C) * delta_t) Javier Civera book P130. A.12
	// q = (q0, q1, q2, q3);
	// output: ( q0 -q1 -q2 -q3
	//			 q1  q0  q3 -q2
	//			 q2 -q3  q0  q1
	//			 q3  q2 -q1  q0)

	double R, X, Y, Z;
	MatrixXd dq3_by_dq2RES(4,4);
	R = q1_in(0);
	X = q1_in(1);
	Y = q1_in(2);
	Z = q1_in(3);
	dq3_by_dq2RES << R, -X, -Y, -Z,
		X,  R,  Z, -Y,
		Y, -Z,  R,  X,
		Z,  Y, -X,  R;
	return dq3_by_dq2RES;
}

MatrixXd MotionModel::dq3_by_dq1(VectorXd q2_in)
{
	// This function calculates dq_k+1^WC/dq((w_k^C + omega^C) * delta_t) which is also 
	// represented by dq3/dq1, see Javier Civera book P132. A.14
	// output: ( q0 -q1 -q2 -q3
	//			 q1  q0 -q3  q2
	//			 q2  q3  q0 -q1
	//			 q3 -q2  q1  q0)
	double R, X, Y, Z;
	MatrixXd dq3_by_dq2RES(4,4);
	R = q2_in(0);
	X = q2_in(1);
	Y = q2_in(2);
	Z = q2_in(3);
	dq3_by_dq2RES << R, -X, -Y, -Z,
		X,  R, -Z,  Y,
		Y,  Z,  R, -X,
		Z, -Y,  X,  R;
	return dq3_by_dq2RES;
}

MatrixXd MotionModel::dqomegadt_by_domega(VectorXd omega, double delta_t)
{
	// This function calculates dq((w_k^C + omega^C) * delta_t)/domega_k^C. Javier book P130. A.17
	// omega Modulus
	double omega_norm = 0;
	size_t omega_size = omega.size(), i;
	MatrixXd dqomegadt_by_domegaRES(4,3);
	for (i = 0; i < omega_size; ++i)
	{
		omega_norm += omega(i) * omega(i);
	}
	omega_norm = sqrt(omega_norm);

	// Use generic ancillary functions to calculate components of Jacobian
	dqomegadt_by_domegaRES << dq0_by_domegaA(omega(0), omega_norm, delta_t), dq0_by_domegaA(omega(1), omega_norm, delta_t), dq0_by_domegaA(omega(2), omega_norm, delta_t),
		dqA_by_domegaA(omega(0), omega_norm, delta_t), dqA_by_domegaB(omega(0), omega(1), omega_norm, delta_t), dqA_by_domegaB(omega(0), omega(2), omega_norm, delta_t),
		dqA_by_domegaB(omega(1), omega(0), omega_norm, delta_t), dqA_by_domegaA(omega(1), omega_norm, delta_t), dqA_by_domegaB(omega(1), omega(2), omega_norm, delta_t), 
		dqA_by_domegaB(omega(2), omega(0), omega_norm, delta_t), dqA_by_domegaB(omega(2), omega(1), omega_norm, delta_t), dqA_by_domegaA(omega(2), omega_norm, delta_t);
	return dqomegadt_by_domegaRES;
}

double MotionModel::dq0_by_domegaA(double omegaA, double omega_norm, double delta_t)
{
	return (-delta_t / 2.0) * (omegaA / omega_norm) * sin(omega_norm * delta_t / 2.0);
}

double MotionModel::dqA_by_domegaA(double omegaA, double omega_norm, double delta_t)
{
	return (delta_t / 2.0) * omegaA * omegaA / (omega_norm * omega_norm) * cos(omega_norm * delta_t / 2.0) 
		+ (1.0 / omega_norm) * (1.0 - omegaA * omegaA / (omega_norm * omega_norm)) * sin(omega_norm * delta_t / 2.0);
}

double MotionModel::dqA_by_domegaB(double omegaA, double omegaB, double omega_norm, double delta_t)
{
	return (omegaA * omegaB / (omega_norm * omega_norm)) * ( (delta_t / 2.0) * cos(omega_norm * delta_t / 2.0)
		- (1.0 / omega_norm) * sin(omega_norm * delta_t / 2.0) );
}




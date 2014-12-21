#include "MathUtil.h"

Vector3d compute_m(double theta, double phi)
{
	double cphi = cos(phi);
	Vector3d m;
	m << cphi*sin(theta), -sin(phi), cphi*cos(theta);
	return m;
}

MatrixXd Q2R(VectorXd Q)
{
	double x, y, z, r;
	MatrixXd R(3,3);
	r = Q(0);
	x = Q(1);
	y = Q(2);
	z = Q(3);
	R << r * r + x * x - y * y - z * z,  2 * (x * y - r * z),            2 * (z * x + r * y),
		 2 * (x * y + r * z),            r * r - x * x + y * y - z * z,  2 * (y * z - r * x),
		 2 * (z * x - r * y),            2 * (y * z + r * x),            r * r - x * x - y * y + z * z;
	return R;
}

MatrixXd dRq_times_a_by_dq(VectorXd q, VectorXd aMat)
{

	MatrixXd dRq_times_a_by_dqRES(3,4), TempR(3,3);
	VectorXd Temp31(3);
	dRq_times_a_by_dqRES.setZero();

	TempR = dR_by_dq0(q);
	Temp31 = TempR * aMat;
	dRq_times_a_by_dqRES.block(0,0,3,1) = Temp31;

	TempR = dR_by_dqx(q);
	Temp31 = TempR * aMat;
	dRq_times_a_by_dqRES.block(0,1,3,1) = Temp31;

	TempR = dR_by_dqy(q);
	Temp31 = TempR * aMat;
	dRq_times_a_by_dqRES.block(0,2,3,1) = Temp31;

	TempR = dR_by_dqz(q);
	Temp31 = TempR * aMat;
	dRq_times_a_by_dqRES.block(0,3,3,1) = Temp31;


	return dRq_times_a_by_dqRES;
}

MatrixXd dR_by_dq0(VectorXd q)
{
	double q0, qx, qy, qz;
	MatrixXd dR_by_dq0RES(3,3);
	q0 = q(0);
	qx = q(1);
	qy = q(2);
	qz = q(3);
	dR_by_dq0RES << 2 * q0, -2 * qz,  2 * qy,
		2 * qz,  2 * q0, -2 * qx,
		-2 * qy,  2 * qx,  2 * q0;
	return dR_by_dq0RES;
}

MatrixXd dR_by_dqx(VectorXd q)
{
	double q0, qx, qy, qz;
	MatrixXd dR_by_dqxRES(3,3);
	q0 = q(0);
	qx = q(1);
	qy = q(2);
	qz = q(3);
	dR_by_dqxRES << 2 * qx,  2 * qy,  2 * qz,
		2 * qy, -2 * qx, -2 * q0,
		2 * qz,  2 * q0, -2 * qx;
	return dR_by_dqxRES;
}

MatrixXd dR_by_dqy(VectorXd q)
{
	double q0, qx, qy, qz;
	MatrixXd dR_by_dqyRES(3,3);
	q0 = q(0);
	qx = q(1);
	qy = q(2);
	qz = q(3);
	dR_by_dqyRES << -2 * qy, 2 * qx,  2 * q0,
		2 * qx, 2 * qy,  2 * qz,
		-2 * q0, 2 * qz, -2 * qy;
	return dR_by_dqyRES;
}

MatrixXd dR_by_dqz(VectorXd q)
{
	double q0, qx, qy, qz;
	MatrixXd dR_by_dqzRES(3,3);
	q0 = q(0);
	qx = q(1);
	qy = q(2);
	qz = q(3);
	dR_by_dqzRES << -2 * qz, -2 * q0, 2 * qx,
		2 * q0, -2 * qz, 2 * qy,
		2 * qx,  2 * qy, 2 * qz;
	return dR_by_dqzRES;
}

// Get points for the ellipse uncertain area
void GetEllPoints(MatrixXd C, VectorXd nu, double chi2, MatrixXd *y)
{
  
  VectorXd th = VectorXd::LinSpaced(Sequential, 100, 0, 2*PI);
  int npoints = th.size();
  VectorXd eig(2,1);
  MatrixXd x_(th.size(), 2), x;
  x_.col(0) = th.array().cos();
  x_.col(1) = th.array().sin();
  x = x_.transpose() * sqrt(chi2);
  EigenSolver<MatrixXd> es(C);
  eig(0) = es.eigenvalues()[0].real();
  eig(1) = es.eigenvalues()[1].real();
  if (eig.minCoeff() <= 0)
    return;
  else
  {
    Eigen::MatrixXd L = C.llt().matrixL();
    MatrixXd K = L.transpose();
    MatrixXd Y(2,npoints);
    Y.row(0) = VectorXd::Ones(npoints).transpose() * nu(0);
    Y.row(1) = VectorXd::Ones(npoints).transpose() * nu(1);
    *y = K * x + Y;
  }
}
#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <math.h>
#include <vector>
#include <string>
#include <Eigen/Eigen>
#include <iostream>

#define PI 3.14159265359
#define eps 2.220446049250313e-016

using namespace std;
using namespace Eigen;

Vector3d compute_m(double theta, double phi);
MatrixXd Q2R(VectorXd Q);

MatrixXd dRq_times_a_by_dq(VectorXd q, VectorXd aMat);
MatrixXd dR_by_dq0(VectorXd q);
MatrixXd dR_by_dqx(VectorXd q);
MatrixXd dR_by_dqy(VectorXd q);
MatrixXd dR_by_dqz(VectorXd q);
void GetEllPoints(MatrixXd C, VectorXd nu, double chi2, MatrixXd *y);

#endif
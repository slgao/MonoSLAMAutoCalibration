#include "Kalman.h"

Kalman::Kalman()
{


}

Kalman::~Kalman()
{


}


// The function predict_camera_measurements and subfunctions are also defined in MapManagementFilterBank file.
// ------------------------------------------------------------------------------------------------------------
void Kalman::predict_camera_measurements(VectorXd x_k_k, KalmanFilter *kalmanfilter, MonoSLAM *mono_slam, Matrix <double, Dynamic, 2> *predicted_measurements)
{

	// This function calculates predicted measurements for each kalman filter
	int size_x_k_k = x_k_k.size(), pos = 0;
	size_t map_features_number = mono_slam->features_info.size(), i;
	// Get the number of predicted measurements for each filter
	kalmanfilter->num_measurements = map_features_number;
	VectorXd cali_para, t_wc, features, yi3d, yi6d, zi;
  // Store reprojection from scopis calibration
  VectorXd zi_scopis;

	MatrixXd r_wc;
	//mono_slam->predicted_measurements = -mono_slam->predicted_measurements.setOnes(map_features_number, 2);
	*predicted_measurements = -(*predicted_measurements).setOnes(map_features_number, 2);
  mono_slam->error_coordinates = -MatrixXd::Ones(map_features_number, 2);

	cali_para = x_k_k.head(kalmanfilter->num_calibration);
	t_wc = x_k_k.segment(kalmanfilter->num_calibration, kalmanfilter->num_t_wc);
	r_wc = Q2R(x_k_k.segment(kalmanfilter->num_calibration + kalmanfilter->num_t_wc, kalmanfilter->num_q_wc));
	features = x_k_k.tail(size_x_k_k - kalmanfilter->num_cam_para);

	for (i = 0; i < map_features_number; ++i)
	{
		if (mono_slam->features_info.at(i).type.compare("cartesian") == 0)
		{
			yi3d = features.segment(pos, 3);
			pos += 3;
			if (hi_cartesian(&zi, yi3d, t_wc, r_wc, mono_slam->Cam, mono_slam->features_info, cali_para))
				(*predicted_measurements).block(i, 0, 1, 2) = zi.transpose();
      if (hi_cartesian(&zi, yi3d, t_wc, r_wc, mono_slam->Cam, mono_slam->features_info, cali_para) && 
        hi_cartesian_scopis(&zi_scopis, yi3d, t_wc, r_wc, mono_slam->Cam, mono_slam->features_info, cali_para))
        mono_slam->error_coordinates.block(i, 0, 1, 2) = zi_scopis.transpose() - zi.transpose();
		}
		if (mono_slam->features_info.at(i).type.compare("inversedepth") == 0)
		{
			yi6d = features.segment(pos, 6);
			pos += 6;
			if (hi_inverse_depth(&zi, yi6d, t_wc, r_wc, mono_slam->Cam, mono_slam->features_info, cali_para))
				(*predicted_measurements).block(i, 0, 1, 2) = zi.transpose();
      if (hi_inverse_depth(&zi, yi6d, t_wc, r_wc, mono_slam->Cam, mono_slam->features_info, cali_para) && 
        hi_inverse_depth_scopis(&zi_scopis, yi6d, t_wc, r_wc, mono_slam->Cam, mono_slam->features_info, cali_para))
        mono_slam->error_coordinates.block(i, 0, 1, 2) = zi_scopis.transpose() - zi.transpose();
		}
	}
}
// ------------------------------------------------------------------------------------------------------------

bool Kalman::hi_cartesian(VectorXd *zi, VectorXd yi3d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
	vector<MonoSLAM::feature_info> features_info, VectorXd cali_para)
{
	VectorXd hrl, uv_u;
	Vector2d uv_d;
	MatrixXd r_cw;
	// Points 3D in camera coordinates
	r_cw = r_wc.inverse();
	hrl = r_cw * (yi3d - t_wc);
	// Check if it is in front of the camera
	if (hrl(2) < 0)
		return false;
	// Project yi3d(camera reference) to image
	uv_u = hu(hrl, cam, cali_para);
	// Add distortion
	uv_d = distort_fm(uv_u, cam, cali_para);

	// Check if the point is in the range of the image
	if ((uv_d(0) > 0 ) && ( uv_d(0) < cam->nCols_ ) && ( uv_d(1)>0 ) && ( uv_d(1) < cam->nRows_))
	{
		* zi = uv_d;
		return true;
	}
	else
		return false;

}

bool Kalman::hi_inverse_depth(VectorXd *zi, VectorXd yi6d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
	vector<MonoSLAM::feature_info> features_info, VectorXd cali_para)
{
	MatrixXd r_cw;
	VectorXd yi, m(3), hrl, uv_u, uv_d;
	double theta, phi, rho;

	r_cw = r_wc.inverse();
	yi = yi6d.head(3);
	theta = yi6d(3);
	phi = yi6d(4);
	rho = yi6d(5);

	m << cos(phi) * sin(theta), -sin(phi), cos(phi) * cos(theta);
	hrl = r_cw * ((yi - t_wc) * rho + m);

	// Check if it is in front of the camera

	if (hrl(2) < 0)
		return false;
	// Project yi6d(camera reference) to image
	uv_u = hu(hrl, cam, cali_para);
	// Add distortion
	uv_d = distort_fm(uv_u, cam, cali_para);
	// Check if the point is in the range of the image
	if (( uv_d(0) > 0 ) && ( uv_d(0) < cam->nCols_ ) && ( uv_d(1) > 0 ) && ( uv_d(1) < cam->nRows_ ))
	{
		*zi = uv_d;
		return true;
	}
	else 
		return false;

}

// Camera reprojection using the calibration parameter got from the scopis calibration method
bool Kalman::hi_cartesian_scopis(VectorXd *zi, VectorXd yi3d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
	vector<MonoSLAM::feature_info> features_info, VectorXd cali_para)
{
	VectorXd hrl, uv_u;
	Vector2d uv_d;
	MatrixXd r_cw;
  // Change calibration parameters to scopis resutls
  cali_para[0] = 440;
  cali_para[1] = 340; //340
  cali_para[2] = 277; //277
  cali_para[3] = 0.0101;
  cali_para[4] = 0.00146;
	// Points 3D in camera coordinates
	r_cw = r_wc.inverse();
	hrl = r_cw * (yi3d - t_wc);
	// Check if it is in front of the camera
	if (hrl(2) < 0)
		return false;
	// Project yi3d(camera reference) to image
	uv_u = hu(hrl, cam, cali_para);
	// Add distortion
	uv_d = distort_fm(uv_u, cam, cali_para);

	// Check if the point is in the range of the image
	if ((uv_d(0) > 0 ) && ( uv_d(0) < cam->nCols_ ) && ( uv_d(1)>0 ) && ( uv_d(1) < cam->nRows_))
	{
		* zi = uv_d;
		return true;
	}
	else
		return false;

}

bool Kalman::hi_inverse_depth_scopis(VectorXd *zi, VectorXd yi6d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
	vector<MonoSLAM::feature_info> features_info, VectorXd cali_para)
{
	MatrixXd r_cw;
	VectorXd yi, m(3), hrl, uv_u, uv_d;
	double theta, phi, rho;
  // Change calibration parameters to scopis resutls
  cali_para[0] = 440;
  cali_para[1] = 340; //340
  cali_para[2] = 277; //277
  cali_para[3] = 0.0101;
  cali_para[4] = 0.00146;

	r_cw = r_wc.inverse();
	yi = yi6d.head(3);
	theta = yi6d(3);
	phi = yi6d(4);
	rho = yi6d(5);

	m << cos(phi) * sin(theta), -sin(phi), cos(phi) * cos(theta);
	hrl = r_cw * ((yi - t_wc) * rho + m);

	// Check if it is in front of the camera

	if (hrl(2) < 0)
		return false;
	// Project yi6d(camera reference) to image
	uv_u = hu(hrl, cam, cali_para);
	// Add distortion
	uv_d = distort_fm(uv_u, cam, cali_para);
	// Check if the point is in the range of the image
	if (( uv_d(0) > 0 ) && ( uv_d(0) < cam->nCols_ ) && ( uv_d(1) > 0 ) && ( uv_d(1) < cam->nRows_ ))
	{
		*zi = uv_d;
		return true;
	}
	else 
		return false;

}
// --------------------------------------------------------------------------------------------------------

VectorXd Kalman::hu(VectorXd yi, Camera *cam, VectorXd cali_para)
{
	double u0, v0, f;
	VectorXd uv_u(2,1);
	u0 = cali_para(1);
	v0 = cali_para(2);
	f = cali_para(0);
	//ku = 1/cam->dx_;
	//kv = 1/cam->dy_;

	uv_u << u0 + (yi(0)/yi(2))*f, v0 + (yi(1)/yi(2))*f;
	return uv_u;
}

Vector2d Kalman::distort_fm(VectorXd uv, Camera *cam, VectorXd cali_para)
{
	double Cx, Cy, k1, k2, dx, dy, xu, yu, ru, rd, f, f_p, D, xd, yd;
	int k;
	Vector2d uvd;

	Cx = cali_para(1);
	Cy = cali_para(2);
	k1 = cali_para(3);
	k2 = cali_para(4);
	dx = cam->dx_;
	dy = cam->dy_;

	xu = (uv(0) - Cx) * dx;
	yu = (uv(1) - Cy) * dy;

	ru = sqrt(xu * xu + yu * yu);
	rd = ru / (1 + k1 * pow(ru, 2) + k2 * pow(ru, 4));

	for (k = 0; k < 25; ++k)
	{
		f = rd + k1 * pow(rd, 3) + k2 * pow(rd, 5) - ru;
		f_p = 1 + 3 * k1 * pow(rd, 2) + 5 * k2 * pow(rd, 4);
		rd = rd - f / f_p;
	}

	D = 1 + k1 * pow(rd, 2) + k2 * pow(rd, 4);
	xd = xu / D;
	yd = yu / D;

	uvd << xd / dx + Cx, yd / dy + Cy;
	return uvd;

}
// ------------------------------------------------------------------------------------------------------------

void Kalman::CalculateDerivatives(VectorXd x_k_km1, MatrixXd PredictedMeasurements, Camera *cam, vector<MonoSLAM::feature_info> features_info, MatrixXd *H_predicted)
{
	// This function calculates the H_predicted, see Javier book P115. 6.4
	size_t num_feat = features_info.size(), // How many features are there in the state vector
		size_x_k_km1 = x_k_km1.size(), num_cam_para = 18, i, pos = 0; // pos is used to select the appropriate element in vector x_features
	VectorXd x_v, x_features, y;
	*H_predicted = MatrixXd::Zero(2 * num_feat, size_x_k_km1); // H_prediceted initialization
	x_v = x_k_km1.head(num_cam_para); // Camera parameters
	x_features = x_k_km1.tail(size_x_k_km1 - num_cam_para); // Feature points
	for (i = 0; i < num_feat; ++i)
	{
		if (PredictedMeasurements(i, 0) != -1.0) // Only work on the matched feature points
		{
			if (features_info.at(i).type.compare("cartesian") == 0)
			{
				y = x_features.segment(pos, 3); // Get the next cartesian feature point
				pos += 3;
				// Calculates Hi derivative cartesian
				CalculateHiCartesian(x_v, y, PredictedMeasurements.row(i), cam, i, features_info, H_predicted);
			}
			else 
			{
				y = x_features.segment(pos, 6); // Get the next inverse depth feature point
				pos += 6;
			    // Calculates Hi derivative inverse depth
				CalculateHiInverseDepth(x_v, y, PredictedMeasurements.row(i), cam, i, features_info, H_predicted);
			}

		}
		else 
		{
			if (features_info.at(i).type.compare("cartesian") == 0)
				pos += 3;
			else
				pos += 6;
		}
	
	}
}

//-----------------------------------------------------------------------------------------------------
///////////////////////                                                /////////////////////////////////
/////////////////////// Calculation of Hi derivative in inverse depth /////////////////////////////////
///////////////////////                                                /////////////////////////////////
//-----------------------------------------------------------------------------------------------------
void Kalman::CalculateHiInverseDepth(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, size_t FeatIndex, 
	vector<MonoSLAM::feature_info> features_info, MatrixXd *H_predicted)
{
	// Computes full Jacobian H, see Javier book P132. A.28 or P141. A.84
	VectorXd calibration = Xv_km1_k.head(5); // Calibration intrinsic parameters
	VectorXi InverseDepthFeatIndex, CartesianFeatIndex;
	cv::KeyPoint feature; // cv::KeyPoint type for the predicted measurements
	feature.pt.x = (float)zi(0); feature.pt.y = (float)zi(1);
	size_t num_feat = features_info.size(), i; 
	InverseDepthFeatIndex = VectorXi::Zero(num_feat);
	CartesianFeatIndex = VectorXi::Zero(num_feat);

	// Here the initialization can be commented out
	// Classify cartesian feature points and inverse depth feature points
	for (i = 0; i < num_feat; ++i)
	{
		if (features_info.at(i).type.compare("inversedepth") == 0)
			InverseDepthFeatIndex(i) = 1;
		else if (features_info.at(i).type.compare("cartesian") == 0)
			CartesianFeatIndex(i) = 1;
	}

	// Initialize H_predicted matrix size with already known cartesian and inverse depth feature points
	(*H_predicted).middleRows(2 * FeatIndex, 2) = MatrixXd::Zero(2, 18 + 3 * CartesianFeatIndex.sum() + 6 * InverseDepthFeatIndex.sum()); 
	// Here the initialization can be commented out

	// Assign H_predicted
	dhd_df(Xv_km1_k, yi, zi, cam, calibration, H_predicted, FeatIndex);
	//cout << "dhd_df" << (*H_predicted).maxCoeff() <<endl;
	dhd_dCx(H_predicted, FeatIndex);
	dhd_dCy(H_predicted, FeatIndex);
	dhd_dk1(cam, feature, calibration, H_predicted, FeatIndex);
	dhd_dk2(cam, feature, calibration, H_predicted, FeatIndex);
	dh_dxv(cam, Xv_km1_k, yi, zi, H_predicted, FeatIndex);

	size_t index_insertion = 18 + 3 * (CartesianFeatIndex.head(FeatIndex)).sum() + 6 * (InverseDepthFeatIndex.head(FeatIndex).sum());
	dh_dy(cam, Xv_km1_k, yi, zi, H_predicted, index_insertion, FeatIndex);
}

void Kalman::dhd_df(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex)
{
	// Computes dhd/df, see Javier book P142. A.92 
	(*H_predicted).block(2 * FeatIndex, 0, 2, 1) = dhd_dhu(cam, zi, calibration) * dhu_df(cam, Xv_km1_k, yi);

}

MatrixXd Kalman::dhd_dhu(Camera *cam, VectorXd zi_d, VectorXd calibration)
{
	// Computes dhd/dhu, see Javier book P133. A.31
	MatrixXd J_undistort(2,2);
	// Calculates the wohle Jacobian
	JacobUndistortFM(cam, zi_d, calibration, &J_undistort);
	return J_undistort.inverse();
}

void Kalman::JacobUndistortFM(Camera *cam, VectorXd uvd, VectorXd calibration, MatrixXd *J_undistort)
{
	// Computes the Jacobian of the undistortion of the image coordinates. See Javier book P133. A.32
	// Input: cam          - camera parameters
	//        calibration  - calibration parameters
	//        uvd          - distorted image points in pixels
	// Output: J_undistort - Jacobian
	double Cx, Cy, k1, k2, dx, dy, ud, vd, xd, yd, rd2, rd4, uu_ud, vu_vd, 
		uu_vd, vu_ud;
	// Extract calibration parameters
	Cx = calibration(1);
	Cy = calibration(2);
	k1 = calibration(3);
	k2 = calibration(4);
	// Get pixel size in x and y direction
	dx = cam->dx_;
	dy = cam->dy_;
	
	ud = uvd(0);
	vd = uvd(1);
	xd = (ud - Cx) * dx;
	yd = (vd - Cy) * dy;

	rd2 = xd * xd + yd * yd;
	rd4 = rd2 * rd2;

	// Calculation of the four elements in the Jacobian matrix
	// uu_ud and vu_vd calculation here is the same as the matlab code and the book

	uu_ud = (1 + k1 * rd2 + k2 * rd4) + (ud - Cx) * (k1 + 2 * k2 * rd2) * (2 * (ud - Cx) * dx * dx);
	vu_vd = (1 + k1 * rd2 + k2 * rd4) + (vd - Cy) * (k1 + 2 * k2 * rd2) * (2 * (vd - Cy) * dy * dy);

	uu_vd = (ud - Cx) * (k1 + 2 * k2 * rd2) * (2 * (vd - Cy) * dy * dy);
	vu_ud = (vd - Cy) * (k1 + 2 * k2 * rd2) * (2 * (ud - Cx) * dx * dx);

	*J_undistort << uu_ud, uu_vd, vu_ud, vu_vd; // 2x2 matrix
}

VectorXd Kalman::dhu_df(Camera *cam, VectorXd Xv_km1_k, VectorXd yi)
{
	// Computes dhu/df, see Javier book P142. A.94. Computation of hu A.23
	double theta, phi, rho, hcx, hcy, hcz;
	VectorXd rw, mi(3), hc, DhuDf(2);
	MatrixXd Rrw;
	//double ku = 1 / cam->dx_;
	//double kv = 1 / cam->dy_;

	rw = Xv_km1_k.segment(5,3);
	// From quaternion to rotation matrix
	Rrw = Q2R(Xv_km1_k.segment(8,4)).inverse();

	theta = yi(3);
	phi = yi(4);
	rho = yi(5);
	mi = compute_m(theta, phi);

	// Features from world reference to camera reference
	hc = Rrw * ((yi.head(3) - rw) * rho + mi); // See Javier book P131. A.21
	hcx = hc(0);
	hcy = hc(1);
	hcz = hc(2);

	//// The expression is the same as the book
	//DhuDf << -hcx / hcz * ku, -hcy / hcz * kv; // 2x1 vector
	// The expression is the same as the matlab code
	DhuDf << hcx / hcz, hcy / hcz; // 2x1 vector
	return DhuDf;
}

void Kalman::dhd_dCx(MatrixXd *H_predicted, size_t FeatIndex)
{
	// Computes dhd/dCx, see Javier book P142. A.95
	// Matlab uses fsolve to calculate J
	VectorXd DhdDCx(2);
	(*H_predicted).block(2 * FeatIndex, 1, 2, 1) << 1.0, 0.0;
}

void Kalman::dhd_dCy(MatrixXd *H_predicted, size_t FeatIndex)
{
	// Computes dhd/dCx, see Javier book P142. A.95
	// Matlab uses fsolve to calculate J 
	VectorXd DhdDCy(2);
	(*H_predicted).block(2 * FeatIndex, 2, 2, 1) << 0.0, 1.0;
}

void Kalman::dhd_dk1(Camera *cam, cv::KeyPoint feature, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex)
{
	// This function is not consistent with the theory explanation
	// in the Javier book P142. A.101
	double U0, V0, k1, k2, ud, vd, xd, yd, rd, uu, vu;
	cv::KeyPoint undistort_feature;
	U0 = calibration(1);
	V0 = calibration(2);
	k1 = calibration(3);
	k2 = calibration(4);

	ud = (double)feature.pt.x;
	vd = (double)feature.pt.y;
	xd = (ud - U0) * cam->dx_;
	yd = (vd - V0) * cam->dy_;
	rd = sqrt(xd * xd + yd * yd);

	// Get the undistorted feature point
	undistort_fm(&undistort_feature, feature, cam, calibration);
	uu = (double)undistort_feature.pt.x;
	vu = (double)undistort_feature.pt.y;
	(*H_predicted).block(2 * FeatIndex, 3, 2, 1) << -(uu * pow(rd, 2)) / pow((1 + k1 * pow(rd, 2) + k2 * pow(rd, 4)), 2), 
		-(vu * pow(rd, 2)) / pow((1 + k1 * pow(rd, 2) + k2 * pow(rd, 4)), 2);
}

void Kalman::dhd_dk2(Camera *cam, cv::KeyPoint feature, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex)
{
	// This function is not consistent with the theory explanation
	// in the Javier book P142. A.101
	double U0, V0, k1, k2, ud, vd, xd, yd, rd, uu, vu;
	cv::KeyPoint undistort_feature;
	U0 = calibration(1);
	V0 = calibration(2);
	k1 = calibration(3);
	k2 = calibration(4);

	ud = (double)feature.pt.x;
	vd = (double)feature.pt.y;
	xd = (ud - U0) * cam->dx_;
	yd = (vd - V0) * cam->dy_;
	rd = sqrt(xd * xd + yd * yd);

	// Get the undistorted feature point
	undistort_fm(&undistort_feature, feature, cam, calibration);
	uu = (double)undistort_feature.pt.x;
	vu = (double)undistort_feature.pt.y;
	(*H_predicted).block(2 * FeatIndex, 4, 2, 1) << -(uu * pow(rd, 4)) / pow((1 + k1 * pow(rd, 2) + k2 * pow(rd, 4)), 2), 
		-(vu * pow(rd, 4)) / pow((1 + k1 * pow(rd, 2) + k2 * pow(rd, 4)), 2);
}

// This undistort_fm method has also been defined in the MapManagementFilterBank.cpp file
// ----------------------------------------------------------------------------------------------------------------
void Kalman::undistort_fm(cv::KeyPoint *undistort_feature, cv::KeyPoint feature, Camera *cam, VectorXd calibration)
{
	double Cx, Cy, k1, k2, dx, dy, xd, yd, rd, D, xu, yu, feature_x, feature_y;
	Cx = calibration(1); Cy = calibration(2); k1 = calibration(3); k2 = calibration(4);
	feature_x = (double)feature.pt.x;
	feature_y = (double)feature.pt.y;
	dx = cam->dx_; dy = cam->dy_;
	xd = (feature_x - Cx) * dx;
	yd = (feature_y - Cy) * dy;
	rd = sqrt(xd * xd + yd * yd);
	D = 1 + k1 * pow(rd, 2) + k2 * pow(rd, 4);
	xu = xd * D;
	yu = yd * D;
	(*undistort_feature).pt.x = (float)(xu/dx + Cx);
	(*undistort_feature).pt.y = (float)(yu/dy + Cy);
}
// ----------------------------------------------------------------------------------------------------------------

void Kalman::dh_dxv(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t FeatIndex)
{
	// This function computes the derivative of each measurement by the camera
	// states except the intrinsic parameters
	// Further the derivative can be seperated into derivative wrt. rW and wrt. Qw.
	// See Javier book P133 A.30
	(*H_predicted).block(2 * FeatIndex, 5, 2, 13) << dh_drw(cam, Xv_km1_k, yi, zi), dh_dqwr(cam, Xv_km1_k, yi, zi), MatrixXd::Zero(2,6);

}


// Calculation of dh_drw
// -------------------------------------------------------------------------------------------------------
MatrixXd Kalman::dh_drw(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi)
{
	// Calculates measurement derivative by rW, see Javier book P133. A.31
	MatrixXd Hi11;
	Hi11 = dh_dhrl(cam, Xv_km1_k, yi, zi) * dhrl_drw(Xv_km1_k, yi);
	return Hi11; // 2x3 matrix
}

MatrixXd Kalman::dh_dhrl(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi)
{
	// dhd/dhu * dhu/dhC
	VectorXd calibration;
	calibration = Xv_km1_k.head(5); // Five intrinsic parameters
	return dhd_dhu(cam, zi, calibration) * dhu_dhrl(cam, Xv_km1_k, yi);
}

MatrixXd Kalman::dhu_dhrl(Camera *cam, VectorXd Xv_km1_k, VectorXd yi)
{
	// This function calculates the derivative for the pinhole camera model, see Javier book P133. A.34
	// Here the expression is not the same as the one in the book
	double f, theta, phi, rho, hcx, hcy, hcz;
	VectorXd rw, mi(3), hc;
	MatrixXd Rrw, a(2,3);

	//double ku = 1 / cam->dx_;
	//double kv = 1 / cam->dy_;
	f = Xv_km1_k(0);
	rw = Xv_km1_k.segment(5,3);
	Rrw = Q2R(Xv_km1_k.segment(8,4)).inverse();

	theta = yi(3);
	phi = yi(4);
	rho = yi(5);
	mi << cos(phi)*sin(theta), -sin(phi), cos(phi)*cos(theta);
	
	// Calculates measurements in camera reference using A.21
	hc = Rrw * ((yi.head(3) - rw) * rho + mi);
	hcx = hc(0);
	hcy = hc(1);
	hcz = hc(2);
	// The form below is the expression in the book
	//a << -f * ku / hcz, 0, hcx * f * ku / pow(hcz, 2), 0, -f * kv / hcz, hcy * f * kv / pow(hcz, 2);
	a << f / hcz, 0, -hcx * f / pow(hcz, 2), 0, f/hcz, -hcy * f / pow(hcz,2);
	return a;
}

MatrixXd Kalman::dhrl_drw(VectorXd Xv_km1_k, VectorXd yi)
{
	// Computes dhC/drW, see Javier book P133. A.35
	return -(Q2R(Xv_km1_k.segment(8,4)).inverse()) * yi(5);
}
// -------------------------------------------------------------------------------------------------------
// Calculation of dh_drw


// Calculation of dh_dqwr
// -------------------------------------------------------------------------------------------------------
MatrixXd Kalman::dh_dqwr(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi)
{
	// Calculates measurement derivative by qW, see Javier book P134. A.37
	MatrixXd Hi12;
	Hi12 = dh_dhrl(cam, Xv_km1_k, yi, zi) * dhrl_dqwr(Xv_km1_k, yi);
	return Hi12; // 2x4 matrix
}

MatrixXd Kalman::dhrl_dqwr(VectorXd Xv_km1_k, VectorXd yi)
{
	// Calculates dhC/dqWC, this can be divided into dhC/dqCW * dqCW/dqWC
	// see Javier book P134. A.38
	double lambda, phi, theta; // lambda is rho
	VectorXd rw, qwr, mi(3), Diag(4);
	rw = Xv_km1_k.segment(5,3);
	qwr = Xv_km1_k.segment(8,4);
	lambda = yi(5);
	phi = yi(4);
	theta = yi(3);
	mi << cos(phi) * sin(theta), -sin(phi), cos(phi) * cos(theta);
	Diag << 1, -1, -1, -1;
	// QConj function is used to transfer qwr from WC to CW, camera reference
	return dRq_times_a_by_dq(QConj(qwr), (yi.head(3) - rw) * lambda + mi) * Diag.asDiagonal();
}

VectorXd Kalman::QConj(VectorXd q)
{
	// Convert the quaternion from world reference to camera reference
	VectorXd q_bar;
	q_bar = -q;
	q_bar(0) = q(0);
	return q_bar;
}
// -------------------------------------------------------------------------------------------------------
// Calculation of dh_dqwr


// Calculation of dh_dy
// -------------------------------------------------------------------------------------------------------
void Kalman::dh_dy(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t index_insertion, size_t FeatIndex)
{
	// Computes the partial derivative dh/dy, it can be divided into dhd/dhu * dhu/dhC * hC/dhy
	// See Javier book P135. A.51
	(*H_predicted).block(2 * FeatIndex, index_insertion, 2, 6) = dh_dhrl(cam, Xv_km1_k, yi, zi) * dhrl_dy(Xv_km1_k, yi);
}

MatrixXd Kalman::dhrl_dy(VectorXd Xv_km1_k, VectorXd yi)
{
	// Calculates dhC/dy, see Javier book P135. A.52
	double lambda, phi, theta;
	VectorXd rw, dmi_dthetai, dm_dtheta(3), dmi_dphii, dm_dphi(3);
	MatrixXd Rrw, a(3,6);
	rw = Xv_km1_k.segment(5,3);
	Rrw = Q2R(Xv_km1_k.segment(8,4)).inverse();
	lambda = yi(5);
	phi = yi(4);
	theta = yi(3);
	
	dm_dtheta << cos(phi) * cos(theta), 0, -cos(phi) * sin(theta);
	dm_dphi << -sin(phi)*sin(theta), -cos(phi), -sin(phi)*cos(theta);
	dmi_dthetai = Rrw * dm_dtheta;
	dmi_dphii = Rrw * dm_dphi;

	a << lambda * Rrw, dmi_dthetai, dmi_dphii, Rrw * (yi.head(3) - rw);
	return a;
}
// -------------------------------------------------------------------------------------------------------
// Calculation of dh_dy
//-----------------------------------------------------------------------------------------------------
///////////////////////                                                /////////////////////////////////
/////////////////////// Calculation of Hi derivative in inverse depth /////////////////////////////////
///////////////////////                                                /////////////////////////////////
//-----------------------------------------------------------------------------------------------------

// ******************************************************************************************************** //

//-----------------------------------------------------------------------------------------------------
///////////////////////                                                /////////////////////////////////
///////////////////////   Calculation of Hi derivative in Cartesian   /////////////////////////////////
///////////////////////                                                /////////////////////////////////
//-----------------------------------------------------------------------------------------------------
void Kalman::CalculateHiCartesian(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, size_t FeatIndex, 
	vector<MonoSLAM::feature_info> features_info, MatrixXd *H_predicted)
{
	// Computes full Jacobian H, see Javier book P132. A.28 or P141. A.84
	VectorXd calibration = Xv_km1_k.head(5); // Calibration intrinsic parameters
	VectorXi InverseDepthFeatIndex, CartesianFeatIndex;
	cv::KeyPoint feature; // cv::KeyPoint type for the predicted measurements
	feature.pt.x = (float)zi(0); feature.pt.y = (float)zi(1);
	size_t num_feat = features_info.size(), i; 
	InverseDepthFeatIndex = VectorXi::Zero(num_feat);
	CartesianFeatIndex = VectorXi::Zero(num_feat);

	// Here the initialization can be commented out
	// Classify cartesian feature points and inverse depth feature points
	for (i = 0; i < num_feat; ++i)
	{
		if (features_info.at(i).type.compare("inversedepth") == 0)
			InverseDepthFeatIndex(i) = 1;
		else if (features_info.at(i).type.compare("cartesian") == 0)
			CartesianFeatIndex(i) = 1;
	}

	// Initialize H_predicted matrix size with already known cartesian and inverse depth feature points
	(*H_predicted).middleRows(2 * FeatIndex, 2) = MatrixXd::Zero(2, 18 + 3 * CartesianFeatIndex.sum() + 6 * InverseDepthFeatIndex.sum()); 
	// Here the initialization can be commented out

	// Assign H_predicted
	dhd_df_xyz(Xv_km1_k, yi, zi, cam, calibration, H_predicted, FeatIndex);
	dhd_dCx(H_predicted, FeatIndex);
	dhd_dCy(H_predicted, FeatIndex);
	dhd_dk1(cam, feature, calibration, H_predicted, FeatIndex);
	dhd_dk2(cam, feature, calibration, H_predicted, FeatIndex);
	dh_dxv_xyz(cam, Xv_km1_k, yi, zi, H_predicted, FeatIndex);

	size_t index_insertion = 18 + 3 * (CartesianFeatIndex.head(FeatIndex)).sum() + 6 * (InverseDepthFeatIndex.head(FeatIndex).sum());
	dh_dy_xyz(cam, Xv_km1_k, yi, zi, H_predicted, index_insertion, FeatIndex);
}

void Kalman::dhd_df_xyz(VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, Camera *cam, VectorXd calibration, MatrixXd *H_predicted, size_t FeatIndex)
{
	// Computes dhd/df, see Javier book P142. A.92 
	// Here yi is a 3d vector, for the reason of cartesian representation
	(*H_predicted).block(2 * FeatIndex, 0, 2, 1) = dhd_dhu(cam, zi, calibration) * dhu_df_xyz(cam, Xv_km1_k, yi);

}

VectorXd Kalman::dhu_df_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi)
{
	// Computes dhu/df, see Javier book P142. A.94. Computation of hu A.23
	double hcx, hcy, hcz;
	VectorXd rw, hc, DhuDf(2);
	MatrixXd Rrw;
	//double ku = 1 / cam->dx_;
	//double kv = 1 / cam->dy_;

	rw = Xv_km1_k.segment(5,3);
	// From quaternion to rotation matrix
	Rrw = Q2R(Xv_km1_k.segment(8,4)).inverse();

	// Features from world reference to camera reference
	hc = Rrw * (yi - rw); // See Javier book P131. A.21
	hcx = hc(0);
	hcy = hc(1);
	hcz = hc(2);

	//// The expression is the same as the book
	//DhuDf << -hcx / hcz * ku, -hcy / hcz * kv; // 2x1 vector
	// The expression is the same as the matlab code
	DhuDf << hcx / hcz, hcy / hcz; // 2x1 vector
	return DhuDf;
}

void Kalman::dh_dxv_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t FeatIndex)
{
	// This function computes the derivative of each measurement by the camera
	// states except the intrinsic parameters
	// Further the derivative can be seperated into derivative wrt. rW and wrt. Qw.
	// See Javier book P133 A.30
	// dh/drw and dh/dqwr are different from the inverse depth version, this calculation is not 
	// contained in the matlab code
	(*H_predicted).block(2 * FeatIndex, 5, 2, 13) << dh_drw_xyz(cam, Xv_km1_k, yi, zi), dh_dqwr_xyz(cam, Xv_km1_k, yi, zi), MatrixXd::Zero(2,6);

}

// Calculation of dh_drw_xyz
// -------------------------------------------------------------------------------------------------------
MatrixXd Kalman::dh_drw_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi)
{
	// Calculates measurement derivative by rW, see Javier book P133. A.31
	MatrixXd Hi11;
	Hi11 = dh_dhrl_xyz(cam, Xv_km1_k, yi, zi) * dhrl_drw_xyz(Xv_km1_k);
	return Hi11; // 2x3 matrix
}

MatrixXd Kalman::dh_dhrl_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi)
{
	// dhd/dhu * dhu/dhC
	VectorXd calibration;
	calibration = Xv_km1_k.head(5); // Five intrinsic parameters
	return dhd_dhu(cam, zi, calibration) * dhu_dhrl_xyz(cam, Xv_km1_k, yi);
}

MatrixXd Kalman::dhu_dhrl_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi)
{
	// This function calculates the derivative for the pinhole camera model, see Javier book P133. A.34
	// Here the expression is not the same as the one in the book
	double f, hcx, hcy, hcz;
	VectorXd rw, mi(3), hc;
	MatrixXd Rrw, a(2,3);

	//double ku = 1 / cam->dx_;
	//double kv = 1 / cam->dy_;
	f = Xv_km1_k(0);
	rw = Xv_km1_k.segment(5,3);
	Rrw = Q2R(Xv_km1_k.segment(8,4)).inverse();
	
	// Calculates measurements in camera reference using A.21
	hc = Rrw * (yi.head(3) - rw);
	hcx = hc(0);
	hcy = hc(1);
	hcz = hc(2);
	// The form below is the expression in the book
	//a << -f * ku / hcz, 0, hcx * f * ku / pow(hcz, 2), 0, -f * kv / hcz, hcy * f * kv / pow(hcz, 2);
	a << f / hcz, 0, -hcx * f / pow(hcz, 2), 0, f/hcz, -hcy * f / pow(hcz,2);
	return a;
}

MatrixXd Kalman::dhrl_drw_xyz(VectorXd Xv_km1_k)
{
	// Computes dhC/drW, see Javier book P133. A.35
	// This is different from the inverse depth version
	return -(Q2R(Xv_km1_k.segment(8,4)).inverse());
}

MatrixXd Kalman::dh_dqwr_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi)
{
	// Calculates measurement derivative by qW, see Javier book P134. A.37
	MatrixXd Hi12;
	Hi12 = dh_dhrl_xyz(cam, Xv_km1_k, yi, zi) * dhrl_dqwr_xyz(Xv_km1_k, yi);
	return Hi12; // 2x4 matrix
}

MatrixXd Kalman::dhrl_dqwr_xyz(VectorXd Xv_km1_k, VectorXd yi)
{
	// Calculates dhC/dqWC, this can be divided into dhC/dqCW * dqCW/dqWC
	// see Javier book P134. A.38
	VectorXd rw, qwr, Diag(4);
	rw = Xv_km1_k.segment(5,3);
	qwr = Xv_km1_k.segment(8,4);
	Diag << 1, -1, -1, -1;
	// QConj function is used to transfer qwr from WC to CW, camera reference
	// Here is different from the inverse depth version
	return dRq_times_a_by_dq(QConj(qwr), (yi.head(3) - rw)) * Diag.asDiagonal();
}
// -------------------------------------------------------------------------------------------------------
// Calculation of dh_drw_xyz

// Calculation of dh_dy_xyz
// -------------------------------------------------------------------------------------------------------
void Kalman::dh_dy_xyz(Camera *cam, VectorXd Xv_km1_k, VectorXd yi, VectorXd zi, MatrixXd *H_predicted, size_t index_insertion, size_t FeatIndex)
{
	// Computes the partial derivative dh/dy, it can be divided into dhd/dhu * dhu/dhC * hC/dhy
	// See Javier book P135. A.51
	(*H_predicted).block(2 * FeatIndex, index_insertion, 2, 3) = dh_dhrl_xyz(cam, Xv_km1_k, yi, zi) * dhrl_dy_xyz(Xv_km1_k);
}

MatrixXd Kalman::dhrl_dy_xyz(VectorXd Xv_km1_k)
{
	// Calculates dhC_xyz/dy, see Javier book P135. A.55
	MatrixXd Rrw;
	Rrw = Q2R(Xv_km1_k.segment(8,4)).inverse();
	return Rrw;
}
// -------------------------------------------------------------------------------------------------------
// Calculation of dh_dy_xyz

//-----------------------------------------------------------------------------------------------------
///////////////////////                                                /////////////////////////////////
///////////////////////   Calculation of Hi derivative in Cartesian   /////////////////////////////////
///////////////////////                                                /////////////////////////////////
//-----------------------------------------------------------------------------------------------------

void Kalman::EKF_Prediction(MonoSLAM *mono_slam)
{
	vector<KalmanFilter *>::iterator FilterBank_iter;
	vector<KalmanFilter *> &FilterBank = mono_slam->KalmanFilterBank_->FilterBank_;
	// Calculation for each filter in the filterbank
	//#pragma omp parallel for
	for (FilterBank_iter = FilterBank.begin(); FilterBank_iter != FilterBank.end(); ++FilterBank_iter)
	{
		mono_slam->MotionModel_->predict_state_and_covariance((*FilterBank_iter)->x_k_k_, (*FilterBank_iter)->p_k_k_, 
			(*FilterBank_iter)->type, (*FilterBank_iter)->std_a_, (*FilterBank_iter)->std_alpha_, &(*FilterBank_iter)->x_k_km1_, &(*FilterBank_iter)->p_k_km1_);

		//(*FilterBank_iter)->predicted_measurements = predict_camera_measurements((*FilterBank_iter)->x_k_km1_, *FilterBank_iter, mono_slam);
		predict_camera_measurements((*FilterBank_iter)->x_k_km1_, *FilterBank_iter, mono_slam, &(*FilterBank_iter)->predicted_measurements);

		CalculateDerivatives((*FilterBank_iter)->x_k_km1_, (*FilterBank_iter)->predicted_measurements, mono_slam->Cam, mono_slam->features_info, &(*FilterBank_iter)->H_predicted_);

		size_t rows_H_predicted = (*FilterBank_iter)->H_predicted_.rows();
		(*FilterBank_iter)->R_predicted_ = MatrixXd::Identity(rows_H_predicted, rows_H_predicted) * (*FilterBank_iter)->std_z_;

		// Calculation of S_predicted, see Javier paper (13), to show the 3sigma area
		(*FilterBank_iter)->S_predicted_ = (*FilterBank_iter)->H_predicted_ * (*FilterBank_iter)->p_k_km1_ * (*FilterBank_iter)->H_predicted_.transpose() + (*FilterBank_iter)->R_predicted_ * 0.25;
		//cout <<"the largest value in S_predicted-------" <<(*FilterBank_iter)->S_predicted_.maxCoeff() << endl;
		(*FilterBank_iter)->Rownum_S_predicted = (*FilterBank_iter)->S_predicted_.rows();
		(*FilterBank_iter)->Colnum_S_predicted = (*FilterBank_iter)->S_predicted_.cols();
	}

}

void Kalman::EKF_Update(vector<KalmanFilter *> FilterBank)
{
	vector<KalmanFilter *>::iterator FilterBank_iter;
	//#pragma omp parallel for
	for (FilterBank_iter = FilterBank.begin(); FilterBank_iter != FilterBank.end(); ++FilterBank_iter)
	{
		Update((*FilterBank_iter)->x_k_km1_, (*FilterBank_iter)->p_k_km1_, (*FilterBank_iter)->H_matching_, 
			(*FilterBank_iter)->R_matching_, (*FilterBank_iter)->z_, (*FilterBank_iter)->h_, &(*FilterBank_iter)->x_k_k_, &(*FilterBank_iter)->p_k_k_);
	}

}

void Kalman::Update(VectorXd x_km1_k, MatrixXd p_km1_k, MatrixXd H, MatrixXd R, VectorXd z, VectorXd h, VectorXd *x_k_k, MatrixXd *p_k_k)
{
	MatrixXd K;
	if (z.size() > 0)
	{
	// Kalman filter gain
	MatrixXd S, Jnorm;
	size_t size_p_k_k;
	S = H * p_km1_k * H.transpose() + R;
	K = p_km1_k * H.transpose() * S.inverse();

	// State vector and covariance matrix update
	*x_k_k = x_km1_k + K * (z - h);
	*p_k_k = (MatrixXd::Identity(p_km1_k.rows(), p_km1_k.rows()) - K * H) * p_km1_k;

	// Normalize quaternion
	Jnorm = NormJac((*x_k_k).segment(8, 4));
	size_p_k_k = (*p_k_k).rows();
	// Change the rows and columns we need to chaneg in the covariance matrix
	(*p_k_k).block(0,8,5,4) = (*p_k_k).block(0,8,5,4) * Jnorm.transpose();
	(*p_k_k).block(5,8,3,4) = (*p_k_k).block(5,8,3,4) * Jnorm.transpose();
	(*p_k_k).block(8,0,4,5) = Jnorm * (*p_k_k).block(8,0,4,5);	
	(*p_k_k).block(8,5,4,3) = Jnorm * (*p_k_k).block(8,5,4,3);	
	(*p_k_k).block(8,8,4,4) = Jnorm * (*p_k_k).block(8,8,4,4) * Jnorm.transpose();	
	(*p_k_k).block(8,12,4,size_p_k_k - 13 + 1) = Jnorm * (*p_k_k).block(8,12,4,size_p_k_k - 13 + 1);	
	(*p_k_k).block(12,8,size_p_k_k - 13 + 1,4) = (*p_k_k).block(12,8,size_p_k_k - 13 + 1,4) * Jnorm.transpose();

	(*x_k_k).segment(8,4) = (*x_k_k).segment(8,4) / (*x_k_k).segment(8,4).norm();
	}
	else
	{
		*x_k_k = x_km1_k;
		*p_k_k = p_km1_k;
		//K << 0;
	}

}

MatrixXd Kalman::NormJac(VectorXd q)
{
	double r, x, y, z;
	MatrixXd J(4,4);
	r = q(0);
	x = q(1);
	y = q(2);
	z = q(3);
	
	J <<  x * x + y * y + z * z,  -r * x,                   -r * y,                   -r * z,
	     -x * r,                   r * r + y * y + z * z,   -x * y,                   -x * z,
	     -y * r,                  -y * x,                    r * r + x * x + z * z,   -y * z,
	     -z * r,                  -z * x,                   -z * y,                    r * r + x * x + y * y;
	J = pow(r * r + x * x + y * y + z * z, -3/2) * J;
	return J;
}

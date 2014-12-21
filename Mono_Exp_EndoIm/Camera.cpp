#include "Camera.h"

Camera::Camera()
    :dx_init(1),
	dy_init(1),
	nRows_init(1),
	nCols_init(1)
{

}
Camera::~Camera()
{

}
void Camera::SetCameraParameters(double dx, double dy, int nRows, int nCols, string model)
{
	dx_ = dx;
	dy_ = dy;
	nRows_ = nRows;
	nCols_ = nCols;
	model_ = model;
	if (dx && dy == 0){
		dx_init = 0;
		dy_init = 0;
	}
}

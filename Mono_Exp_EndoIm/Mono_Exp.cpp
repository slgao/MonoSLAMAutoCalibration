// Monoslam_experiment.cpp : Defines the entry point for the console application.

#include "Monoslam.h"
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <vector>
#include <iomanip>

using namespace cv;

int main(int argc, char* argv[])
{
	FILE *Pfile;
	double dx, dy, error_type_I = 0.05, error_type_II = 0.1; //error_type_I = 0.01, error_type_II = 0.05; 
	int nRows, nCols; 
	int frame_id = 0; Frame frame_last, frame_next;
	int initIm = 1, lastIm = 900, step, step_ = 344; // Variable initIm starts at 1 good results step_ = 210, 344, 64, 85
	string model; string CamDataPath = "camdata/cam.txt";

	MonoSLAM *mono_slam = NULL;
	mono_slam = new MonoSLAM();
	mono_slam->Init(CamDataPath, &dx, &dy, &nRows, &nCols, &model, initIm, lastIm);

	mono_slam->Init_filterbank(mono_slam->KalmanFilterBank_);
	//mono_slam->frame_grabber_->GetFrame(frame_id, &frame_last);
	//mono_slam->frame_grabber_->GetFrame(frame_id + 1, &frame_next);
	// ---------------------------------------------------------------------
	// Read image manually
	char s[200]; // Here s size is limited
	//sprintf(s, "images/rawoutput%04d.pgm", step_);
  sprintf_s(s, 200, "exp_scopis_im_2_720x576/step%d.jpg", step_);
  //sprintf_s(s, 200, "images2/seqtest%04d.ppm", step_);
	frame_last.data = imread(s, 0);

  //// Add mask to the image
  //Frame maskedIm_last, maskedIm_next;
  //Frame Im_resize_last, Im_resize_next;
  //cv::Mat mask(frame_last.data.size(), frame_last.data.type());
  //mask.setTo(cv::Scalar(0,0,0));  

  //cv::Point center(160, 120);
  //int radius = 150;
  //cv::circle(mask, center, radius, cv::Scalar(5, 5, 5), -1, 8, 0);
  //frame_last.data.copyTo(maskedIm_last.data,mask);

  //cv::resize(maskedIm_last.data, Im_resize_last.data, cv::Size2i(maskedIm_last.data.cols/2, maskedIm_last.data.rows/2));
 //

  
	sprintf_s(s, 200, "exp_scopis_im_2_720x576/step%d.jpg", step_ + 1);
  //sprintf_s(s, 200, "images2/seqtest%04d.ppm", step_ + 1);
	frame_next.data = imread(s, 0);

  //frame_next.data.copyTo(maskedIm_next.data,mask);
  //cv::resize(maskedIm_next.data, Im_resize_next.data, cv::Size2i(maskedIm_next.data.cols/2, maskedIm_next.data.rows/2));



  if (!frame_last.data.data || !frame_next.data.data)
  {
    return -1;
  }
  
  //imshow("iamge", frame_next.data);
  //cv::waitKey(0);
	// ---------------------------------------------------------------------
	//camera parameter initialization
	// Set camera calibration
	mono_slam->Cam->SetCameraParameters(dx, dy, nRows, nCols, model);
	mono_slam->likelihood_ratio = VectorXd::Zero(mono_slam->KalmanFilterBank_->filter_size);
	mono_slam->A = log10((1 - error_type_II) / error_type_I);
    mono_slam->B = log10(error_type_II / (1 - error_type_I));
	//cout << "dx is " << mono_slam->Cam->dx_ << "\ndy is " << mono_slam->Cam->dy_ <<"\nnRows is "
	//	<< mono_slam->Cam->nRows_ << "\nnCols is " << mono_slam->Cam->nCols_ << "\nmodel is " <<
	//	mono_slam->Cam->model_ << "\ninput_mode is " << mono_slam->input_mode << endl;
	//mono_slam->KalmanFilterBank_->initialize_filterbank();

	fopen_s(&Pfile, "cali_results.txt","a");
	for (step = step_; step < lastIm; ++step)
	{

		mono_slam->GoOneStep(step, initIm, lastIm, frame_last, frame_next);

    //maskedIm_last.data = maskedIm_next.data.clone();
    frame_next.data.copyTo(frame_last.data);
		//frame_last.data = frame_next.data;
		//mono_slam->frame_grabber_->GetFrame(step + 1, &frame_next);
		// ------------------------------------------------------------------
		sprintf_s(s, 200, "exp_scopis_im_2_720x576/step%d.jpg", step + 2);
    //sprintf_s(s, 200, "images2/seqtest%04d.ppm", step + 2);
		frame_next.data = imread(s, 0);

    //frame_next.data.copyTo(maskedIm_next.data,mask);
    //cv::resize(maskedIm_next.data, Im_resize_next.data, cv::Size2i(maskedIm_next.data.cols/2, maskedIm_next.data.rows/2));


    //imshow("iamge", frame_next.data);
    //cv::waitKey(0);
		// ------------------------------------------------------------------
		//fprintf(Pfile, "Step %d: -----------------------------\n %lf\n %lf\n %lf\n %lf\n %lf\n %d\n", step, mono_slam->x_k_k_output(0), mono_slam->x_k_k_output(1),
    //	mono_slam->x_k_k_output(2),mono_slam->x_k_k_output(3),mono_slam->x_k_k_output(4), mono_slam->KalmanFilterBank_->filter_size);

    // Make sure there are always features in the map
    if (mono_slam->features_info.size() == 0)
    {
      continue;
      cv::imshow("image 1", frame_next.data);
    }

    double std_f = std::sqrt(mono_slam->p_k_k_output(0,0));
    double std_Cx = std::sqrt(mono_slam->p_k_k_output(1,1));
    double std_Cy = std::sqrt(mono_slam->p_k_k_output(2,2));
    double std_k1 = std::sqrt(mono_slam->p_k_k_output(3,3));
    double std_k2 = std::sqrt(mono_slam->p_k_k_output(4,4));
    double f = mono_slam->x_k_k_output(0);
    double Cx = mono_slam->x_k_k_output(1);
    double Cy = mono_slam->x_k_k_output(2);
    double k1 = mono_slam->x_k_k_output(3);
    double k2 = mono_slam->x_k_k_output(4);

    fprintf(Pfile, "%d           %lf           %lf            %lf            %lf            %lf            %d            %f            %f            %f            %f            %f            %f            %f            %f            %f            %f\n", 
      step, mono_slam->x_k_k_output(0), mono_slam->x_k_k_output(1),
      mono_slam->x_k_k_output(2),mono_slam->x_k_k_output(3),mono_slam->x_k_k_output(4), mono_slam->KalmanFilterBank_->filter_size, f - 3* std_f, f + 3* std_f, Cx - 3 * std_Cx, Cx + 3 * std_Cx,
      Cy - 3 * std_Cy, Cy + 3 * std_Cy, k1 - 3 * std_k1, k1 + 3 * std_k1, k2 - 3 * std_k2, k2 + 3 * std_k2);
		//printf("Step %d: -----------------------------\n %lf\n %lf\n %lf\n %lf\n %lf\n %d\n", step, mono_slam->x_k_k_output(0), mono_slam->x_k_k_output(1),
		//	mono_slam->x_k_k_output(2),mono_slam->x_k_k_output(3),mono_slam->x_k_k_output(4), mono_slam->KalmanFilterBank_->filter_size);
		fclose(Pfile);
		fopen_s(&Pfile, "cali_results.txt","a");

    
    cout << setprecision(8) <<"--------------------------" << endl <<step << endl 
      << "Focal Length:" <<mono_slam->x_k_k_output(0) << "    " << 194 - 3 * std_f << "    " << 194 + 3 * std_f << endl 
      << "Cx:          " <<mono_slam->x_k_k_output(1) << "    " << 160.2 - 3 * std_Cx << "    " << 160.2 + 3 * std_Cx << endl 
      << "Cy:          " <<mono_slam->x_k_k_output(2) << "    " << 128.9 - 3 * std_Cy << "    " << 128.9 + 3 * std_Cy << endl 
      << setprecision(6)
      << "Kappa1:      " << mono_slam->x_k_k_output(3) * pow(0.0112,2)<< "    " << 0.0633 - 3 * std_k1 << "    " << 0.0633 + 3 * std_k1 << endl 
      << "Kappa2:      " << mono_slam->x_k_k_output(4) * pow(0.0112,4)<< "    " << 0.0139 - 3 * std_k2 << "    " << 0.0139 + 3 * std_k2 << endl 
      << "Filter Size: " <<mono_slam->KalmanFilterBank_->filter_size <<endl;
	}
  
	//cout << "the No. of the frame is " << mono_slam->initialized_features.front().index_frame <<endl << 
	//	"the feature points are " <<endl<< mono_slam->initialized_features.front().initialized_features_per_frame <<endl;

    //namedWindow( "Display window", WINDOW_AUTOSIZE );// Create a window for display.
	//imshow( "Display window", frame.data.rowRange(100,200).colRange(100,300));

	/*
	// -----------------------------------------  OpenCV_FAST ALGORITHM ------------------------------------ //
	std::vector<KeyPoint> keyPoints;
	//cv::FASTX(frame.data, keyPoints, 50, true, 2);
	FastFeatureDetector fast(50); // define detector threshold, TYPE_5_8 = 0, TYPE_7_12 = 1, TYPE_9_16 = 2(default)
    // feature point detection
    fast.detect(frame.data,keyPoints);
    drawKeypoints(frame.data, keyPoints, frame.data, Scalar::all(255), DrawMatchesFlags::DRAW_OVER_OUTIMG);
    imshow("FAST feature", frame.data);
	// -----------------------------------------  OpenCV_FAST ALGORITHM ------------------------------------ //
	*/
  delete mono_slam;
	cv::waitKey(0);
	//cout <<"the size of the kalman filter is "<< KalmanFilterBank -> FilterBank_.size() << endl;
	//cout << KalmanFilterBank->FilterBank_[10]->H_matching_ << endl;
	cout << "END OF THE PROGRAMME!!" <<endl;
	system("pause>nul");
	return 0;
}


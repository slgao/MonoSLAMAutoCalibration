//#include "Monoslam.h"
// Important! without it the MapManagement can not be recognized, while there is no MapManagement header file included in the Monoslam.h
#include "MapManagementFilterBank.h"
#include "Kalman.h"
#include "MotionModel.h"
#include "DataAssociator.h"
#include <time.h>  

MonoSLAM::MonoSLAM():min_number_of_features_in_image(85)
{
	KalmanFilterBank_ = NULL;
	Cam = NULL;
	MapManager = NULL;
	Kalman_ = NULL;
	MotionModel_ = NULL;
	DataAssociator_ = NULL;
}
MonoSLAM::~MonoSLAM()
{
	if (KalmanFilterBank_ != NULL)
		delete  KalmanFilterBank_;

	if (Cam != NULL)
		delete  Cam;

	if (MapManager != NULL)
		delete  MapManager;

	if (Kalman_ != NULL)
		delete  Kalman_;

	if (MotionModel_ != NULL)
		delete  MotionModel_;

	if (DataAssociator_ != NULL)
		delete  DataAssociator_;

  while (!initialized_features.empty())
    initialized_features.pop_back();

  while (!features_info.empty())
    features_info.pop_back();
}
void MonoSLAM::Init(const string &config_path, double* dx, double* dy,
	int* nRows, int* nCols, string* model, int init_frame_id, int init_max)
{

	KalmanFilterBank_ = new FilterBank();
	Cam = new Camera();

	//frame_grabber_ = new FrameGrabber(); // Commented out to release mode
	ParseVarsFile(config_path, dx, dy, nRows, nCols, model, &input_mode, input_name);
	// get frame range from init_frame_id to init_max
	//frame_grabber_->Init(input_name, input_mode, init_frame_id, init_max); // Commented out to release mode 
	MapManager = new MapManagement();
	Kalman_ = new Kalman();
	MotionModel_ = new MotionModel();
	DataAssociator_ = new DataAssociator();
}

void MonoSLAM::Init_filterbank(FilterBank *KalmanFilterBank_)
{
	KalmanFilterBank_->initialize_x_and_p();
	x_k_k_output = KalmanFilterBank_->x_k_k;
	p_k_k_output = KalmanFilterBank_->p_k_k;
	KalmanFilterBank_->initialize_filterbank();
	mu = VectorXd::Ones(KalmanFilterBank_->filter_size) / KalmanFilterBank_->filter_size;
}

void MonoSLAM::GoOneStep(int step, int init_frame_id, int init_max, Frame frame_last, Frame frame_next)
{
	double normalization_coefficient = 0, normalization_coefficient_data = 0;
	VectorXd nu;
	MatrixXd mu_data, S;
	size_t num_filters;
  double start1 = clock();
	MapManager->map_management_filter_bank(step, frame_last, this);
  double end1 = clock();
  // If there is no feature in the map after processing the image, jump out
  if (this->features_info.size() == 0)
  {
    return;
  }
	//for (size_t i = 0; i < features_info.size(); ++i)
	//{
	//	cout << features_info.at(i).uv_when_initialized <<endl;
	//}
	//frame_grabber_->GetFrame(step + 1, &frame_last);

  double start2 = clock();
  Kalman_->EKF_Prediction(this);
  double end2 = clock();
  double start3 = clock();
	DataAssociator_->FilterBankMatching(this, this->mu, frame_next);
  double end3 = clock();
  double start4 = clock();
	Kalman_->EKF_Update(this->KalmanFilterBank_->FilterBank_);
  double end4 = clock();


  cout << double((end1-start1)) / CLOCKS_PER_SEC << "    " << double((end2-start2)) / CLOCKS_PER_SEC << "    " 
    << double((end3-start3)) / CLOCKS_PER_SEC << "    "  << double((end4-start4)) / CLOCKS_PER_SEC << "    " << endl;
	// Compute a posteriori probability for each filter
	num_filters = this->KalmanFilterBank_->filter_size ;
	mu_data = VectorXd::Zero(num_filters);

	for (size_t i = 0; i < num_filters; ++i)
	{
		nu = this->KalmanFilterBank_->FilterBank_.at(i)->z_ - this->KalmanFilterBank_->FilterBank_.at(i)->h_; // Javier Book P114 6.2 innovations vector
		S = this->KalmanFilterBank_->FilterBank_.at(i)->S_matching_; // Javier Book P114 6.2 Covariance of innovations vector
		double s = sqrt(pow(2 * PI, S.rows()) * S.determinant());
		double e = -0.5 * nu.transpose() * S.inverse() * nu;
		double ee = std::exp(e);
		mu_data(i) = 1.0 / s * ee;
		this->mu(i) *= mu_data(i);
		normalization_coefficient_data += mu_data(i);
		normalization_coefficient += this->mu(i);
	}

  mu_data = mu_data / normalization_coefficient_data;
  mu = mu / normalization_coefficient;

	// Final output of state vector and covariance matrix
	this->x_k_k_output = VectorXd::Zero(this->KalmanFilterBank_->FilterBank_.at(0)->x_k_k_.size());
	this->p_k_k_output = MatrixXd::Zero(this->KalmanFilterBank_->FilterBank_.at(0)->p_k_k_.rows(), this->KalmanFilterBank_->FilterBank_.at(0)->p_k_k_.cols());

	// State vector computation
	for (size_t i = 0; i < num_filters; ++i)
	{
		this->x_k_k_output += this->KalmanFilterBank_->FilterBank_.at(i)->x_k_k_ * this->mu(i);
	}

	// Covariance matrix compuatation
	for (size_t i = 0; i < num_filters; ++i)
	{
		this->p_k_k_output += (this->KalmanFilterBank_->FilterBank_.at(i)->p_k_k_ + 
			(this->x_k_k_output - this->KalmanFilterBank_->FilterBank_.at(i)->x_k_k_) * 
			((this->x_k_k_output - this->KalmanFilterBank_->FilterBank_.at(i)->x_k_k_).transpose())) * this->mu(i);
	}


	// Do hypothesis test to see which filters are to be pruned

	if (num_filters > 1)
	{
		VectorXi filters_to_prune;
		filters_to_prune = VectorXi::Zero(num_filters);
		for (size_t i = 0; i < num_filters; ++i)
		{
			this->likelihood_ratio(i) = this->likelihood_ratio(i) + log10(mu_data(i) / ((mu_data.sum() - mu_data(i)) / (num_filters - 1)));
			if (this->likelihood_ratio(i) < B)
				filters_to_prune(i) = 1;
		}

		// Prune bad filters
		size_t num_filters_new;
		num_filters_new = num_filters - filters_to_prune.sum();
		VectorXd mu_new(num_filters_new), likelihood_ratio_new(num_filters_new);
		size_t count = 0, shift = 0;
		for (size_t i = 0; i < num_filters; ++i)
		{
			if (filters_to_prune(i) == 1)
			{
				size_t pos = i - shift;
				this->KalmanFilterBank_->FilterBank_.erase(this->KalmanFilterBank_->FilterBank_.begin() + pos);
				shift ++;
			}
			else
			{
				mu_new(count) = this->mu(i);
				likelihood_ratio_new(count) = this->likelihood_ratio(i);
				count ++;
			}
		}

		// Assign new mu and likelihood ratio
    this->mu = mu_new;
    this->likelihood_ratio = likelihood_ratio_new;
    this->KalmanFilterBank_->filter_size = num_filters_new;

    // Renormalize probabilities
    double sum_of_discrete_probabilities = this->mu.sum();
    this->mu /= sum_of_discrete_probabilities;

  }
  
  // Display the matched and unmatched feature points
  cv::Mat ImPre, ImNPre, Imfinal;
  cv::KeyPoint feature_P, feature_NP, feature_Meas;
  vector<cv::KeyPoint> P_pre, P_npre, P_Meas;
  static int fontFace = cv::FONT_HERSHEY_SCRIPT_SIMPLEX;
  static double fontScale = 0.3;
  static int thickness = 1;
  static struct feat_indices
  {
    string num;
    cv::Point position;
    bool predicted;
  } feat_in;
  vector<feat_indices> ind_feats;
  static struct uncertain
  {
    cv::Scalar color;
    MatrixXd EllPoints;
  } Uncer;
  vector<uncertain> Predicted, Unpredicted;

  int pre_size = this->predicted_measurements.rows();
  for (int index = 0; index < pre_size; ++index)
  {
    if (this->measurements(index, 0) != -1 && measurements(index, 1) != -1)
    {
      feature_P.pt.x = this->predicted_measurements(index, 0);
      feature_P.pt.y = this->predicted_measurements(index, 1);
      feature_P.size = 3;
      P_pre.push_back(feature_P);
      feature_Meas.pt.x = measurements(index, 0);
      feature_Meas.pt.y = measurements(index, 1);
      feature_Meas.size = 3;
      P_Meas.push_back(feature_Meas);
      feat_in.num = to_string((_ULonglong)this->features_info.at(index).ind_feat);
      feat_in.position = cv::Point(feature_Meas.pt.x + 5, feature_Meas.pt.y - 2);
      feat_in.predicted = true;
      ind_feats.push_back(feat_in);

      Uncer.color = cv::Scalar(255, 0, 255);
      MatrixXd Ellpoint;
      double chi2 = 9.2103;
      GetEllPoints(this->S_predicted_.block(2 * index, 2 * index, 2, 2), 
        predicted_measurements.row(index), chi2, &(Uncer.EllPoints));
      Predicted.push_back(Uncer);
    }
    else
    {
      feature_NP.pt.x = this->predicted_measurements(index, 0);
      feature_NP.pt.y = this->predicted_measurements(index, 1);
      feature_NP.size = 3;
      P_npre.push_back(feature_NP);
      feat_in.num = to_string((_ULonglong)this->features_info.at(index).ind_feat);
      feat_in.position = cv::Point(feature_NP.pt.x + 5, feature_NP.pt.y + 2);
      feat_in.predicted = false;
      ind_feats.push_back(feat_in);

      Uncer.color = cv::Scalar(255, 255, 0);
      MatrixXd Ellpoint;
      double chi2 = 9.2103;
      GetEllPoints(this->S_predicted_.block(2 * index, 2 * index, 2, 2), 
        predicted_measurements.row(index), chi2, &(Uncer.EllPoints));
      Unpredicted.push_back(Uncer);
    }
  }
        cv::Mat frame_color;
        cv::cvtColor(frame_next.data, frame_color, CV_GRAY2RGB);
        // draw keyPoints on the image
        for (auto it = Unpredicted.begin(), itEnd = Unpredicted.end(); it != itEnd; it ++)
        {
          for (int j = 0; j < (*it).EllPoints.cols() - 1; j ++)
          {
            cv::Point pt1, pt2;
            pt1 = cv::Point((*it).EllPoints(0,j),(*it).EllPoints(1,j));
            pt2 = cv::Point((*it).EllPoints(0,j+1),(*it).EllPoints(1,j+1));
            cv::line(frame_color, pt1, pt2, CV_RGB(0,0,250), 1, 8);
          }
        }
        for (auto it = Predicted.begin(), itEnd = Predicted.end(); it != itEnd; it ++)
        {
          for (int j = 0; j < (*it).EllPoints.cols() - 1; j ++)
          {
            cv::Point pt1, pt2;
            pt1 = cv::Point((*it).EllPoints(0,j),(*it).EllPoints(1,j));
            pt2 = cv::Point((*it).EllPoints(0,j+1),(*it).EllPoints(1,j+1));
            cv::line(frame_color, pt1, pt2, CV_RGB(250,0,0), 1, 8);
          }
        }

        //ImPre = frame_next.data.clone();
        /*cv::drawKeypoints(frame_next.data, P_pre, ImPre, cv::Scalar(255, 0, 255), cv::DrawMatchesFlags::DEFAULT);
        cv::drawKeypoints(ImPre, P_npre, ImNPre, cv::Scalar(255, 255, 0), cv::DrawMatchesFlags::DEFAULT);*/
        cv::drawKeypoints(frame_color, P_Meas, ImPre, cv::Scalar(0, 0, 255), 4);
        cv::drawKeypoints(ImPre, P_pre, ImNPre, cv::Scalar(0, 255, 0), 4);
        cv::drawKeypoints(ImNPre, P_npre, Imfinal, cv::Scalar(255, 0, 0), 4);
        for (auto it = ind_feats.begin(), itEnd = ind_feats.end(); it != itEnd; it ++)
        {
          if (!(*it).predicted)
          cv::putText(Imfinal, (*it).num, (*it).position, fontFace, fontScale,
          CV_RGB(0,0,250), thickness, 8);
          else
          cv::putText(Imfinal, (*it).num, (*it).position, fontFace, fontScale,
          CV_RGB(250,0,0), thickness, 8);
        }
        //char imname1[200], imname2[200];
        //sprintf(imname1, "%d__1", step);
        //sprintf(imname2, "%d__2", step);
        //imshow(imname1, ImPre);
        //imshow(imname2, ImNPre);
        // Draw endoscopic circle on the image
        int radius = 350; // 330
        cv::Point_<int> center(352, 255); // 360, 288
        const cv::Scalar color(255, 255, 255);
        cv::circle(Imfinal, center, radius, color, 1, 8, 0);
        cv::imshow("image 1", Imfinal);
        
        // Plot error coordinates
        double std_f = std::sqrt(this->p_k_k_output(0,0));

        static cv::vector<cv::KeyPoint> error;
        error.clear();
        static cv::KeyPoint errorP;
        static cv::Mat Imerror(200,200,CV_8UC1);
        Imerror.setTo(cv::Scalar(255, 255, 255));
        errorP.pt.x = 100;
        errorP.pt.y = 100;
        errorP.size = 10;
        error.push_back(errorP);
        FILE *ProjectionError;
        char *errorname = new char[50];
        sprintf(errorname, "errors/error%04d.txt", step);
        fopen_s(&ProjectionError, errorname, "a");
        for (int ind = 0; ind < this->error_coordinates.rows(); ind++)
        {
          if (this->error_coordinates(ind, 0) != -1 && this->error_coordinates(ind, 1) != -1 && this->measurements(ind, 0) != -1 && measurements(ind, 1) != -1)
          {
            fprintf(ProjectionError, "%lf    %lf\n", this->error_coordinates(ind, 0) + 100, this->error_coordinates(ind, 1) + 100);
            errorP.pt.x = this->error_coordinates(ind, 0) + 100;
            errorP.pt.y = this->error_coordinates(ind, 1) + 100;
            errorP.size = 2;
            if (errorP.pt.x > 0 && errorP.pt.x < 200 & errorP.pt.y > 0 && errorP.pt.y < 200)
              error.push_back(errorP);
          }
        }
        fclose(ProjectionError);
        delete[] errorname;
        cv::drawKeypoints(Imerror, error, Imerror, cv::Scalar(0, 0, 100), 4);
        cv::imshow("error", Imerror);
        static char s[30];
        sprintf(s, "error/error%d.jpg", step);
        //cv::imwrite(s, Imerror);


        // Plot pre and measurement error coordinates
        MatrixXd error_p_m;
        error_p_m = -MatrixXd::Ones(this->measurements.rows(), 2);
        for (int ind = 0; ind < this->measurements.rows(); ind++)
        {
          if (this->measurements(ind, 0) != -1 && measurements(ind, 1) != -1)
          {
            error_p_m(ind, 0) = this->predicted_measurements(ind, 0) - this->measurements(ind, 0);
            error_p_m(ind, 1) = this->predicted_measurements(ind, 1) - this->measurements(ind, 1);
          }
        }

        static cv::vector<cv::KeyPoint> errvec_p_m;
        errvec_p_m.clear();
        static cv::KeyPoint errorP_p_m;
        static cv::Mat Imerror_p_m(200,200,CV_8UC1);
        Imerror_p_m.setTo(cv::Scalar(255, 255, 255));
        errorP_p_m.pt.x = 100;
        errorP_p_m.pt.y = 100;
        errorP_p_m.size = 10;
        errvec_p_m.push_back(errorP_p_m);
        for (int ind = 0; ind < error_p_m.rows(); ind++)
        {
          if (error_p_m(ind,0) != -1)
          {
          errorP_p_m.pt.x = error_p_m(ind, 0) + 100;
          errorP_p_m.pt.y = error_p_m(ind, 1) + 100;
          errorP_p_m.size = 2;
          if (errorP_p_m.pt.x > 0 && errorP_p_m.pt.x < 200 & errorP_p_m.pt.y > 0 && errorP_p_m.pt.y < 200)
           errvec_p_m.push_back(errorP_p_m);
          }
        }

        //cv::drawKeypoints(Imerror_p_m, errvec_p_m, Imerror_p_m, cv::Scalar(0, 0, 200), 4);
        //cv::imshow("error_p_m", Imerror_p_m);

        sprintf(s, "error_p_m/error_p_m%d.jpg", step);
        //cv::imwrite(s, Imerror_p_m);


        //static char s[50];
        sprintf(s, "imwrite/step%d.jpg", step);
        cv::imwrite(s, Imfinal);
        cv::waitKey(1);



}


#include "DataAssociator.h"


DataAssociator::DataAssociator()
{


}


DataAssociator::~DataAssociator()
{

}

void DataAssociator::FilterBankMatching(MonoSLAM *mono_slam, VectorXd mu_k_km1, Frame frame)
{
  size_t i, j, k;
  vector<KalmanFilter *> &FilterBank = mono_slam->KalmanFilterBank_->FilterBank_;
	size_t num_filters = mono_slam->KalmanFilterBank_->filter_size;
	size_t num_predicted_measurements = FilterBank.at(0)->num_measurements;
	MatrixXi predicted_by_all_filters, predicted_by_filter(num_predicted_measurements, 2);
	MatrixXd predicted, predicted_transpose, S_predicted, predicted_measurement_transpose; // For the fusing EKF filter bank
	MatrixXd measurements;
	VectorXd z, h; // Vector form of measurements and predicted measurements
	VectorXi matched_index;
	if (num_filters > 1)
	{
		predicted_by_all_filters = MatrixXi::Ones(num_predicted_measurements, 2); // num_pre_mea x 2 one matrix

		// Find the predicted measurements which are successfully predicted by all of the filters
		// The result is a binary matrix of size num_predicted_measurements x 2, 0 value means the point is not 
		// successfully predicted by each of the filters.
		for (i = 0; i < num_filters; ++i)
		{
			// Mark the successfully predicted measurements by 1. Loop over column first because it's column wise stored
			for (k = 0; k < 2; ++k)
			{
				for (j = 0; j < num_predicted_measurements; ++j)
					predicted_by_filter(j,k) = FilterBank.at(i)->predicted_measurements(j,k) > 0 ? 1 : 0;
			}
			// Check whether it is a successful prediction when the value > 0
			for (k = 0; k < 2; ++k)
			{
				for (j = 0; j < num_predicted_measurements; ++j)
					predicted_by_all_filters(j,k) = predicted_by_all_filters(j,k) && predicted_by_filter(j,k);
			}
			
		}

		// Using the mixed EKF filter bank prediction to match
		// Fusing EKF filter bank
		predicted = MatrixXd::Zero(num_predicted_measurements, 2);
		S_predicted = MatrixXd::Zero(FilterBank.at(0)->Rownum_S_predicted, FilterBank.at(0)->Colnum_S_predicted);

		for (i = 0; i < num_filters; ++i)
		{
			// Accumulates the predicted for every filter
			predicted +=  FilterBank.at(i)->predicted_measurements * mu_k_km1(i);
		}

		// Here we use transpose and map it into a vector for the sake of calculating the S_predicted matrix below
		predicted_transpose = predicted.transpose();
		Map<VectorXd> predicted_transpose_vec(predicted_transpose.data(), predicted_transpose.rows() * predicted_transpose.cols());
		for (i = 0; i < num_filters; ++i)
		{
			// Accumulates S_predicted for every filter
			predicted_measurement_transpose = FilterBank.at(i)->predicted_measurements.transpose();
			// Reshape predicted measurements, first transpose because of the column storage
			Map<VectorXd> predicted_measure_transpose_vec(predicted_measurement_transpose.data(), 
				predicted_measurement_transpose.rows() * predicted_measurement_transpose.cols());
			MatrixXd matt = (predicted_transpose_vec - predicted_measure_transpose_vec) * 
				((predicted_transpose_vec - predicted_measure_transpose_vec).transpose());
			S_predicted += (FilterBank.at(i)->S_predicted_ + (predicted_transpose_vec - predicted_measure_transpose_vec) * 
				(predicted_transpose_vec - predicted_measure_transpose_vec).transpose()) * mu_k_km1(i);
		}
		// Remove the unpredicted measurements by all filters
		for (k = 0; k < 2; ++k)
		{
			for (j = 0; j < num_predicted_measurements; ++j)
			{
				if (predicted_by_all_filters(j,k) == 0)
					predicted(j,k) = -1;
			}
		}
	}
	else
	{
		predicted = FilterBank.at(0)->predicted_measurements;
		S_predicted = FilterBank.at(0)->S_predicted_;
	}
	
	Matching(frame, predicted, S_predicted, mono_slam->features_info, mono_slam->Cam, &measurements, &z, &h, &matched_index);

	//JointCompatibility(predicted, measurements, S_predicted); // not implemented yet

	// fill in the filter banks
	size_t size_z = z.size(), count;
	for (i = 0; i < num_filters; ++i)
	{
		FilterBank.at(i)->h_ = h;
		FilterBank.at(i)->z_ = z;
		FilterBank.at(i)->measurements_ = measurements;
		FilterBank.at(i)->R_matching_ = MatrixXd::Identity(size_z, size_z) * FilterBank.at(i)->std_z_;
		FilterBank.at(i)->H_matching_.resize(matched_index.size() * 2, FilterBank.at(i)->H_predicted_.cols());
		count = 0;
		for (j = 0; j < (size_t)matched_index.size(); ++j)
		{
			FilterBank.at(i)->H_matching_.middleRows(count * 2, 2) = FilterBank.at(i)->H_predicted_.middleRows(matched_index(j) * 2, 2);
			count ++;
		}
		FilterBank.at(i)->S_matching_ = FilterBank.at(i)->H_matching_ * FilterBank.at(i)->p_k_km1_ * FilterBank.at(i)->H_matching_.transpose() + 
			MatrixXd::Identity(FilterBank.at(i)->H_matching_.rows(), FilterBank.at(i)->H_matching_.rows()) * FilterBank.at(i)->std_z_;
	}
	mono_slam->measurements = measurements;
	mono_slam->predicted_measurements = predicted;
  mono_slam->S_predicted_ = S_predicted; // Ellipse region
}

void DataAssociator::Matching(Frame frame, MatrixXd predicted_measurements, MatrixXd S_predicted, vector<MonoSLAM::feature_info> features_info, Camera *cam, MatrixXd *measurements,
	 VectorXd *measurements_map, VectorXd *predicted_measurements_map, VectorXi *matched_index)
{
	double correlation_threshold = 0.85, chi2inv_2_95 = 5.99146; // the value of chi1inv_2_95 is obtained directly from the matlab
  int i, j, k, jj, kk, l, m;
	size_t num_features, half_patch_size_when_matching, pixels_in_the_match_patch, index_predicted, num_matched = 0;
	num_features = predicted_measurements.rows();
	half_patch_size_when_matching = features_info.at(0).half_patch_size_when_matching;
	pixels_in_the_match_patch = (size_t)pow((double)2 * half_patch_size_when_matching + 1, 2);
	MatrixXd S, invS;
	(*measurements).resize(num_features, 2);
	(*measurements).fill(-1);

	index_predicted = 0;

	for (i = 0; i < (int)num_features; ++i)
	{
		if (predicted_measurements(i,0) != -1 && predicted_measurements(i,1) != -1)
		{
			index_predicted ++;
			S = S_predicted.block(2 * i, 2 * i, 2, 2);
			if (sqrt(S.diagonal().maxCoeff()) < 100)
			{
				invS = S.inverse();
        cv::Mat predicted_patch;
        features_info.at(i).patch_when_initialized.copyTo(predicted_patch);
				//cv::Mat predicted_patch = features_info.at(i).patch_when_initialized;
				size_t half_search_region_size_x = (size_t)ceil(2 * sqrt(S(0,0)));
				size_t half_search_region_size_y = (size_t)ceil(2 * sqrt(S(1,1)));
        // make sure the half_search_region is not too large
        //if (half_search_region_size_x > 15)
        //{
        //  half_search_region_size_x = 15;
        //}
        //else if (half_search_region_size_y > 15)
        //{
        //  half_search_region_size_y = 15;
        //}

				//Map<MatrixXf> predicted_patch(((features_info.at(i)).patch_when_initialized).data());
				MatrixXd patches_for_correlation = MatrixXd::Zero(pixels_in_the_match_patch, (2 * half_search_region_size_x + 1) * (2 * half_search_region_size_y + 1) + 1);
				MatrixXd match_candidates = MatrixXd::Zero(2,(2 * half_search_region_size_x + 1) * (2 * half_search_region_size_y + 1) + 1);
				predicted_patch.convertTo(predicted_patch, CV_64FC1);
				int corr_row_index = 0; // Row index for each column of patches_for_correlation
				for (j = 0; j < predicted_patch.cols; ++j)
					for (k = 0; k < predicted_patch.rows; ++k)
					{
						patches_for_correlation(corr_row_index,0) = predicted_patch.at<double>(j,k);
						corr_row_index ++;
					}
				size_t index_patches_for_correlation = 0;
        for (jj = (int)(floor(predicted_measurements(i,0) + 0.5) - half_search_region_size_x); jj <  (int)(floor(predicted_measurements(i,0) + 0.5) + half_search_region_size_x + 1); ++jj)
        {
          for (kk = (int)(floor(predicted_measurements(i,1) + 0.5) - half_search_region_size_y); kk <  (int)(floor(predicted_measurements(i,1) + 0.5) + half_search_region_size_y + 1); ++kk)
          {
            VectorXd nu(2);
            nu << jj - predicted_measurements(i,0), kk - predicted_measurements(i,1);
            if ((nu.transpose() * invS * nu) < chi2inv_2_95)
            {
              if ((jj > (int)half_patch_size_when_matching) && (jj < cam->nCols_ - (int)half_patch_size_when_matching) &&
                (kk > (int)half_patch_size_when_matching) && (kk < cam->nRows_ - (int)half_patch_size_when_matching))
              {
                // Copy the subimage to image_patch but point to the position of the ROI
                cv::Mat image_patch;
                cv::Rect rect(jj - half_patch_size_when_matching, kk - half_patch_size_when_matching, 
                  2 * half_patch_size_when_matching + 1, 2 * half_patch_size_when_matching + 1);
                frame.data(rect).copyTo(image_patch);

                //image_patch = frame.data.rowRange(kk - half_patch_size_when_matching - 1, kk + half_patch_size_when_matching).
                //	colRange(jj - half_patch_size_when_matching - 1, jj + half_patch_size_when_matching);
                // rowRange and colRange method is changed to cv:.Rect method
                //image_patch = frame.data(cv::Rect(jj - half_patch_size_when_matching, kk - half_patch_size_when_matching, 
                //  2 * half_patch_size_when_matching + 1, 2 * half_patch_size_when_matching + 1));
                index_patches_for_correlation ++;
                image_patch.convertTo(image_patch, CV_64FC1);
                corr_row_index = 0;
                for (l = 0; l < image_patch.rows; ++l)
                  for (m = 0; m < image_patch.cols; ++m)
                  {
                    patches_for_correlation(corr_row_index,index_patches_for_correlation) = image_patch.at<double>(l,m);
                    corr_row_index ++;
                  }
                  match_candidates(0, index_patches_for_correlation - 1) = jj;
                  match_candidates(1, index_patches_for_correlation - 1) = kk;
              }
            }
          }
        }

        // If no match_candidates there will be no measured points
        if (index_patches_for_correlation != 0)
        {
          int r, c;
          double maxcorr;

          // Calculates the correlation matrix
          MatrixXd patches_for_correlation_block; 
          patches_for_correlation_block = patches_for_correlation.block(0,0,pixels_in_the_match_patch,index_patches_for_correlation + 1);
          VectorXd Mean_correlation;
          Mean_correlation = patches_for_correlation_block.colwise().sum() / pixels_in_the_match_patch; // The mean value for different image patches
          patches_for_correlation_block.rowwise() -= Mean_correlation.transpose(); // Eigen broadcasting, each row subtracts the mean value
          MatrixXd Cov;
          Cov = patches_for_correlation_block.transpose() * patches_for_correlation_block / (pixels_in_the_match_patch - 1);
          VectorXd Cov_diagonal;
          Cov_diagonal = Cov.diagonal();
          //for (size_t ind = 0; ind <  Cov_diagonal.size(); ++ind)
          //	Cov_diagonal(ind) = sqrt(Cov_diagonal(ind));
          Cov_diagonal = Cov_diagonal.cwiseSqrt();
          MatrixXd Cov_diagonal_square;
          Cov_diagonal_square = Cov_diagonal * Cov_diagonal.transpose();
          //Cov = Cov.array() / Cov_diagonal_square.array();

          // We can reduce the computational cost further, e.g. just compute the first row value

          Cov = Cov.cwiseQuotient(Cov_diagonal_square); // Coefficient-wise division gets the correlation matrix

          //for (int index = 0; index < Cov.rows(); ++index)
          //	Cov(index,index)  = 0;
          //cout << Cov.maxCoeff() <<endl << Cov.block(0,1,1,Cov.cols() - 1).maxCoeff() << endl;
          //for (size_t index = 0; index < Cov.cols(); ++index)
          //{
          //	if (_isnan(Cov(0,index)))
          //		Cov(0,index) = 0;
          //}
          maxcorr = Cov.block(0,1,1,Cov.cols() - 1).maxCoeff(&r, &c); // Find the maximum value and row and column index of it

          if (maxcorr > correlation_threshold){
            VectorXd match_can = match_candidates.middleCols(c, 1);
            if (sqrt((double)((match_can(0)-352) * (match_can(0)-352) + (match_can(1)-255) * (match_can(1)-255))) < 350 -10)
            (*measurements).middleRows(i, 1) = match_candidates.middleCols(c, 1).transpose(); // Fill in the matched feature points
          }
        }
      }
    }
  }
  for (int index = 0; index <  (*measurements).rows(); ++index)
  {
    if ((*measurements)(index,0) != -1)
      num_matched ++;
  }

  MatrixXd measurements_update(num_matched, 2), predicted_measurements_update(num_matched, 2);
  (*matched_index).resize(num_matched); // Stores the matched features indices
  int count = 0;
  for (int index = 0; index <  (*measurements).rows(); ++index)
  {
    if ((*measurements)(index,0) != -1)
    {
      measurements_update(count,0) = (*measurements)(index,0);
      measurements_update(count,1) = (*measurements)(index,1);
      predicted_measurements_update(count,0) = predicted_measurements(index,0);
      predicted_measurements_update(count,1) = predicted_measurements(index,1);
      (*matched_index)(count) = index;
      count ++;
    }
  }
  MatrixXd measurements_update_transpose = measurements_update.transpose();
  MatrixXd predicted_measurements_update_transpose = predicted_measurements_update.transpose();
  Map<VectorXd> measurements_map_(measurements_update_transpose.data(), num_matched * 2);
	Map<VectorXd> predicted_measurements_map_(predicted_measurements_update_transpose.data(), num_matched * 2);
	*measurements_map = measurements_map_;
	*predicted_measurements_map = predicted_measurements_map_;
}

void DataAssociator::JointCompatibility(MatrixXd predicted, MatrixXd measurements, MatrixXd S_predicted)
{
	
}

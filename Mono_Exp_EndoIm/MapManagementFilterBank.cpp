#include "MapManagementFilterBank.h"

MapManagement::MapManagement(): max_attempts(25), max_initialization_attempts(25), 
  half_patch_size_when_initialized(18), half_patch_size_when_matching(18), imsize_determined(false), 
  horizontal_div(6), vertical_div(6)
{

}

MapManagement::~MapManagement()
{


}

void MapManagement::update_features_info(MonoSLAM* mono_slam)
{
  int number_of_features;
  number_of_features = mono_slam->measurements.rows();
  for (int i = 0; i < number_of_features; ++i)
  {
    if (mono_slam->predicted_measurements(i,1) != -1)
      mono_slam->features_info.at(i).times_predicted += 1;
    if (mono_slam->measurements(i,1) != -1 && mono_slam->measurements(i,1) != -2)
      mono_slam->features_info.at(i).times_measured += 1;	
  }

}

void MapManagement::inversedepth_2_cartesian_filter_bank(MonoSLAM* mono_slam)
{
  int num_features = mono_slam -> features_info.size();
  if (num_features == 0)
    return;
  int num_filters = mono_slam-> KalmanFilterBank_ ->FilterBank_.size();
  double linearity_index_threshold = 0.1; 
  double std_rho, rho, std_d, theta, phi, d_c2p, cos_alpha, linearity_index; 
  VectorXd X, x_filter_i_tail, x_filter_i_head; 
  MatrixXd P, J(3,6); MatrixXd Id = MatrixXd::Identity(3,3);
  X = mono_slam -> x_k_k_output; P = mono_slam -> p_k_k_output;
  int size_reduce = 3, size_x_filter_i, size_X_old = X.size();
  MatrixXd J_all(size_X_old - 3, size_X_old);
  Vector3d m, x_c1, x_c2, p, dm_dtheta, dm_dphi;
  vector<MonoSLAM::feature_info>& features_info = mono_slam->features_info;
  vector<KalmanFilter*>& filter_bank = mono_slam-> KalmanFilterBank_ ->FilterBank_;
  for (int i = 1; i < num_features; ++i)
  {
    if (features_info.at(i).type == "inversedepth")
    {
      // we assume the basic state vector(without features) has a length of 18
      int initialPositionOfFeature = 18;
      for (int j = 0; j < i; ++j)
      {
        if (features_info.at(j).type == "cartesian") initialPositionOfFeature += 3;
        if (features_info.at(j).type == "inversedepth") initialPositionOfFeature += 6;
      }
      std_rho = sqrt(P(initialPositionOfFeature + 5, initialPositionOfFeature + 5));
      rho = X(initialPositionOfFeature + 5);
      std_d = std_rho/(pow(rho,2));
      theta = X(initialPositionOfFeature + 3);
      phi = X(initialPositionOfFeature + 4);
      m = compute_m(theta, phi);
      x_c1 = X.segment(initialPositionOfFeature, 3);
      x_c2 = X.segment(5,3);
      p = inversedepth2cartesian(X.segment(initialPositionOfFeature, 6));
      Vector3d p_ = p - x_c1, p__ = p - x_c2;
      d_c2p = p__.norm();
      cos_alpha = (double)(p_.transpose() * p__ )/ (p_.norm() * p__.norm());
      linearity_index = 4*std_d*cos_alpha/d_c2p;
      if (linearity_index < linearity_index_threshold)
      {
        for (int i_filters = 0; i_filters < num_filters; ++i_filters)
        {
          VectorXd x_filter_i = filter_bank.at(i_filters) ->x_k_k_;
          MatrixXd p_filter_i = filter_bank.at(i_filters) ->p_k_k_;
          p = inversedepth2cartesian(x_filter_i.segment(initialPositionOfFeature,6));
          theta = x_filter_i(initialPositionOfFeature + 3);
          phi = x_filter_i(initialPositionOfFeature + 4);
          rho = x_filter_i(initialPositionOfFeature + 5);

          // update x_filter_i
          size_x_filter_i = x_filter_i.size();
          x_filter_i_head = x_filter_i.head(initialPositionOfFeature);
          x_filter_i_tail = x_filter_i.tail(size_x_filter_i - initialPositionOfFeature - 6);
          x_filter_i.resize(size_x_filter_i - size_reduce);
          x_filter_i << x_filter_i_head, p, x_filter_i_tail;
          dm_dtheta << cos(phi)*cos(theta), 0, -cos(phi)*sin(theta);
          dm_dphi << -sin(phi)*sin(theta), -cos(phi),  -sin(phi)*cos(theta);
          J << Id, (1/rho) * dm_dtheta, (1/rho) * dm_dphi, -m/(pow(rho, 2));
          J_all << MatrixXd::Identity(initialPositionOfFeature, initialPositionOfFeature), MatrixXd::Zero(initialPositionOfFeature, 6), MatrixXd::Zero(initialPositionOfFeature, size_X_old - initialPositionOfFeature - 6),
            MatrixXd::Zero(3, initialPositionOfFeature), J, MatrixXd::Zero(3, size_X_old - initialPositionOfFeature - 6), 
            MatrixXd::Zero(size_X_old - initialPositionOfFeature - 6, initialPositionOfFeature), MatrixXd::Zero(size_X_old - initialPositionOfFeature - 6, 6),
            MatrixXd::Identity(size_X_old - initialPositionOfFeature - 6, size_X_old - initialPositionOfFeature - 6);
          //p_filter_i.resize(J_all.cols(),J_all.rows());
          p_filter_i = J_all * p_filter_i * J_all.transpose();
          filter_bank.at(i_filters) -> x_k_k_ = x_filter_i;
          filter_bank.at(i_filters) -> p_k_k_ = p_filter_i;
        }
        features_info.at(i).type = 	"cartesian";
        return; // one step one conversion
      }

    }
  }

}

Vector3d MapManagement::inversedepth2cartesian(VectorXd inverse_depth)
{
  Vector3d rw, m, cartesian; double theta, phi, rho, cphi;
  rw = inverse_depth.segment(0,3);
  theta = inverse_depth(3);
  phi = inverse_depth(4);
  rho = inverse_depth(5);

  cphi = cos(phi);
  m << cphi*sin(theta), -sin(phi), cphi*cos(theta);
  cartesian(0) = rw(0) + (1/rho) * m(0);
  cartesian(1) = rw(1) + (1/rho) * m(1);
  cartesian(2) = rw(2) + (1/rho) * m(2);
  return cartesian;
}

void MapManagement::delete_features_filter_bank(MonoSLAM* mono_slam)
{
  int num_filters = mono_slam -> KalmanFilterBank_ ->FilterBank_.size(), i, j, num_features = mono_slam -> features_info.size();
  vector<MonoSLAM::feature_info>& features_info = mono_slam->features_info;
  vector<KalmanFilter*>& filter_bank = mono_slam -> KalmanFilterBank_ ->FilterBank_;
  vector<int> deletion_list;
  for (i = 1; i < num_features; ++i)
  {
    if ((features_info.at(i).times_measured < 0.5 * features_info.at(i).times_predicted) &&
      (features_info.at(i).times_predicted > 5))
      deletion_list.push_back(i);
  }
  if (deletion_list.size() != 0)
  {
    for (i = deletion_list.size() - 1; i >= 0; i= i - 1)
    {
      for (j = 0; j <  num_filters; ++j)
      {
        delete_a_feature(filter_bank.at(j) -> x_k_k_, filter_bank.at(j) -> p_k_k_, deletion_list.at(i), features_info);
      }
      features_info.erase(features_info.begin() + deletion_list.at(i));
    }
  }
}

void MapManagement::delete_a_feature(VectorXd &x_k_k, MatrixXd &p_k_k, int FeatToDelete,vector<MonoSLAM::feature_info> features_info)
{
  int parToDelete, indexToDelete, i;
  if (features_info.at(FeatToDelete).type.compare("cartesian") == 0)
    parToDelete = 3;
  else parToDelete = 6;
  // we assume there are five calibration parameters
  indexToDelete = 18;

  for (i = 0; i < FeatToDelete - 1; ++i)
  {
    if (features_info.at(i).type.compare("inversedepth") == 0)
      indexToDelete += 6;
    else if (features_info.at(i).type.compare("cartesian") == 0)
      indexToDelete += 3;
  }

  // update x_k_k
  int size_x_k_k = x_k_k.size();
  VectorXd x_k_k_head = x_k_k.head(indexToDelete);
  VectorXd x_k_k_tail = x_k_k.tail(size_x_k_k - indexToDelete - parToDelete);
  x_k_k.resize(size_x_k_k - parToDelete);
  x_k_k << x_k_k_head, x_k_k_tail;

  // update p_k_k block wise
  int size_p_k_k = p_k_k.rows();
  int size_p_k_k_tail = size_p_k_k - indexToDelete - parToDelete;
  MatrixXd p_k_k_block00 = p_k_k.block(0,0,indexToDelete,indexToDelete);
  MatrixXd p_k_k_block01 = p_k_k.block(0,indexToDelete + parToDelete, indexToDelete, size_p_k_k_tail);
  MatrixXd p_k_k_block10 = p_k_k.block(indexToDelete + parToDelete, 0, size_p_k_k_tail, indexToDelete);
  MatrixXd p_k_k_block11 = p_k_k.block(indexToDelete + parToDelete, indexToDelete + parToDelete, size_p_k_k_tail, size_p_k_k_tail);
  p_k_k.resize(size_p_k_k - parToDelete, size_p_k_k - parToDelete);
  p_k_k << p_k_k_block00, p_k_k_block01, p_k_k_block10, p_k_k_block11;
}

void MapManagement::map_management_filter_bank(int step, Frame frame, MonoSLAM *mono_slam)
{
  /*int num_filters = mono_slam->KalmanFilterBank_->FilterBank_.size();*/
  int num_of_matched = 0;
  size_t i;
  bool at_least_one_initialized = false;
  // Update features info
  update_features_info(mono_slam);

  // Convert features from inverse depth to cartesian, if necessary
  inversedepth_2_cartesian_filter_bank(mono_slam);

  // Delete features, if necessary
  delete_features_filter_bank(mono_slam);

  // Initialize features randomly
  int measurements_size = mono_slam->measurements.rows();
  if (measurements_size == 0)
    initialize_features_filter_bank(step, mono_slam->min_number_of_features_in_image - measurements_size, mono_slam, frame);
  else
  {
    find_matched_measurements( &num_of_matched, mono_slam->measurements);
    if (num_of_matched < mono_slam->min_number_of_features_in_image)
      initialize_features_filter_bank(step, mono_slam->min_number_of_features_in_image - num_of_matched, mono_slam, frame);
  }

  // Save the feature points
  MonoSLAM::uv_initialized UV_initialized;
  UV_initialized.features_size = 0;
  UV_initialized.index_frame = step;
  for (i = 0; i <  mono_slam->features_info.size(); ++i)
  {
    if (mono_slam->features_info.at(i).init_frame == step)
    {
      UV_initialized.features_size ++;
      at_least_one_initialized = true;
      UV_initialized.isempty = false;
    }
  }
  UV_initialized.initialized_features_per_frame.resize(2, UV_initialized.features_size);
  int feat_per_frame_count = 0;
  for (i = 0; i <  mono_slam->features_info.size(); ++i)
  {
    if (mono_slam->features_info.at(i).init_frame == step)
    {
      UV_initialized.initialized_features_per_frame(0, feat_per_frame_count) = mono_slam->features_info.at(i).init_measurement(0);
      UV_initialized.initialized_features_per_frame(1, feat_per_frame_count) = mono_slam->features_info.at(i).init_measurement(1);
      feat_per_frame_count ++;
    }
  }

  if (at_least_one_initialized == false)
  {
    UV_initialized.isempty = true;
    mono_slam->initialized_features.push_back(UV_initialized);
  }
  else
    mono_slam->initialized_features.push_back(UV_initialized);
}

void MapManagement::find_matched_measurements(int *num_of_matched, MatrixXd measurements)
{   
  // Finds the number of matched feature points
  for (int i = 0; i < measurements.rows(); ++i)
  {
    if (measurements(i,0) != -1.0)
      (*num_of_matched) ++;
  }
}

void MapManagement::initialize_features_filter_bank(int step, int num_feat_to_init, MonoSLAM *mono_slam, Frame frame)
{
  int attempts = 0, initialized = 0;
  cv::KeyPoint feature_TBadded;
  // Pointer reference to mono_slam->kalmanFitlerBank_
  FilterBank *filterbank = mono_slam->KalmanFilterBank_;

  // use just one Filter to calculate the predicted_measurements for the sake of computational reduce
  mono_slam->predicted_measurements = predict_camera_measurements(filterbank->FilterBank_.back()->x_k_k_, filterbank->FilterBank_.back(), mono_slam); 
  while((initialized <  num_feat_to_init) && (attempts < max_attempts))
  {
    attempts ++;
    feature_added = false;
    initialize_a_feature_other_filter_bank(step, frame, mono_slam, &feature_TBadded);
    if (feature_added)
      initialized ++;
    mono_slam->predicted_measurements = predict_camera_measurements(filterbank->FilterBank_.back()->x_k_k_, filterbank->FilterBank_.back(), mono_slam);
  }
}

void MapManagement::initialize_a_feature_other_filter_bank(int step, Frame frame, MonoSLAM *mono_slam, cv::KeyPoint *feature_TBadded)
{

   //Circle detector for endoscopic image 
   //Reduce the noise so we avoid false circle detection
    //cv::Mat blurIm;
    //GaussianBlur( frame.data, blurIm, cv::Size(9, 9), 2, 2 );
    ////declare and initialize both parameters that are subjects to change
    //int cannyThreshold = 110;
    //int accumulatorThreshold = 39;
    //// runs the actual detection
    //std::vector<cv::Vec3f> circles;
    //HoughCircles( blurIm, circles, CV_HOUGH_GRADIENT, 1, frame.data.rows/8, cannyThreshold, accumulatorThreshold, 0, 0 );

    //// Discard unplausible circles
    //if (circles.size() > 1)
    //{
    //  
    //}
    
    //// clone the input image for displaying purposes
    //cv::Mat display = frame.data.clone();
    //for( size_t i = 0; i < circles.size(); i++ )
    //{
    //cv::Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
    //int radius = cvRound(circles[i][2]);
    //// circle center
    //circle( display, center, 3, cv::Scalar(0,255,0), -1, 8, 0 );
    //// circle outline
    //circle( display, center, radius, cv::Scalar(255,255,255), 3, 8, 0 );
    //}

    //// shows the results
    //const std::string windowName = "Hough Circle Detection Demo";
    //cv::namedWindow( windowName, cv::WINDOW_AUTOSIZE );
    //imshow( windowName, display);
    //cv::waitKey(0);
  
  // Reference to FilterBank_
  int i_, i, j, k, l, m, *index_row_min = new int[vertical_div], 
    *row_min = new int[vertical_div], index_col_min, col_min, num_keyPoints, 
    rand_index, index_feat_to_be_added, filter_bank_size;
  int left_bound, right_bound, up_bound, down_bound, hori_interval, verti_interval;
  bool detected_new = false;
  double initial_rho = 1, std_rho = 1;
  VectorXi hori_boundaries(horizontal_div + 1), verti_boundaries(vertical_div + 1);
  VectorXd newFeature(6);
  vector<KalmanFilter*> &filter_bank = mono_slam->KalmanFilterBank_->FilterBank_;
  filter_bank_size = filter_bank.size();
  MatrixXi initialized;
  // uv_pred is the form * *
  //                     * *
  //	    			         * *
  Matrix <double, Dynamic, 2>	&uv_pred = mono_slam->predicted_measurements;
  Matrix <double, Dynamic, 2> NewkeyPoints;
  if (!imsize_determined)
  {
    cv::Size S = frame.data.size();
    horizontal_size = S.width;
    vertical_size = S.height;
  }

  // Set up horizontal boundaries
  left_bound = half_patch_size_when_initialized;
  right_bound = horizontal_size - half_patch_size_when_initialized - 1;
  hori_interval = (right_bound - left_bound) / horizontal_div;
  for (i = 0; i < horizontal_div; ++i)
  {
  hori_boundaries(i) = left_bound + i * hori_interval;
  }
  hori_boundaries(horizontal_div) = right_bound;

  // Set up vertical boundaries
  up_bound = half_patch_size_when_initialized;
  down_bound = vertical_size - half_patch_size_when_initialized - 1;
  verti_interval = (down_bound - up_bound) / vertical_div;
  for (i = 0; i < vertical_div; ++i)
  {
  verti_boundaries(i) = up_bound + i * verti_interval;
  }
  verti_boundaries(vertical_div) = down_bound;

  //////
  //hori_boundaries.resize(4);
  //verti_boundaries.resize(4);
  //hori_boundaries << 6, 105, 211, 312;
  //verti_boundaries << 6, 79, 159, 232;
  //////

  initialized.setZero(horizontal_div, vertical_div);
  if (uv_pred.size() == 0)
    uv_pred = VectorXd::Zero(2).transpose(); // Set uv_pred in row order
  int uv_pred_rows = uv_pred.rows();
  int uv_pred_cols = uv_pred.cols();
  Matrix<bool, Dynamic, Dynamic> uv_pred_flag(uv_pred_rows, 2);

  // Search the image area which has the less initialized points
  for (i = 0; i < horizontal_div; ++i)
  {
    for (j = 0; j < vertical_div; ++j)
    {

      int sum = 0;
      // Loop over uv_pred_flag elements to check points whether in the current image block wrt. x or y coordinates
      for (k = 0; k < uv_pred_cols; ++k)
      {
        for (l = 0; l < uv_pred_rows; ++l)
        {
          uv_pred_flag(l, k) = (k == 0) ? 
            uv_pred(l, k) > hori_boundaries(i) && uv_pred(l, k) < hori_boundaries(i + 1) : 
          uv_pred(l, k) > verti_boundaries(j) && uv_pred(l, k) < verti_boundaries(j + 1);
        }
      }
      for (m = 0; m< uv_pred_rows; ++m)
      {
        // Change in place, 1 at the 0th column means the points locates in the current image block
        uv_pred_flag(m, 0) = uv_pred_flag(m, 0) && uv_pred_flag(m, 1);
        if (uv_pred_flag(m, 0))
          sum ++;
      }
      initialized(j, i) = sum; // Position of the initialized elements represents the position of the image block 
    }
  }
  for (i_ = 0; i_ < max_initialization_attempts; ++i_)
  {
    if (detected_new)
      break;

    // Search the minimum in each row of the initialized matrix
    for (j = 0; j < vertical_div; ++j)
    {
      index_row_min[j] = 0;
      row_min[j] = initialized(j, 0);
      for (k = 0; k < horizontal_div; ++k)
      {
        if (row_min[j] > initialized(j, k))
        {
          row_min[j] = initialized(j, k);
          index_row_min[j] = k;
        }
      }
    }

    index_col_min = 0;
    col_min = row_min[0];
    for (l = 0; l < vertical_div; ++l)
    {
      if (col_min > row_min[l])
      {
        col_min = row_min[l];
        index_col_min = l;
      }
    }

    // index_col_min is the index of row in initialized --- index_row_min[index_col_min] is the index of column in initialized
    // in case 1 1 0(*)
    //		     0 0 0
    //         0 0 0 index_col_min is 0, index_row_min[index_col_min] is 2

    // Get image block where feature points are to be found, here we should copy the subimage but point to the ROI
    cv::Mat im_block;
    cv::Rect rect((int)(hori_boundaries(index_row_min[index_col_min])), 
      (int)(verti_boundaries(index_col_min)), 
      hori_boundaries(index_row_min[index_col_min] + 1) - hori_boundaries(index_row_min[index_col_min]) + 1, 
      verti_boundaries(index_col_min + 1) - verti_boundaries(index_col_min) + 1);
    frame.data(rect).copyTo(im_block);

    //imshow("sd", im_block);
    //cv::waitKey(0);

    //cv::Mat im_block = frame.data.rowRange(verti_boundaries(index_col_min) - 1, verti_boundaries(index_col_min + 1)).
    //	colRange(hori_boundaries(index_row_min[index_col_min]) - 1, hori_boundaries(index_row_min[index_col_min] + 1));

    //cv::Mat im_block = frame.data(cv::Rect((int)(hori_boundaries(index_row_min[index_col_min])), 
    //  (int)(verti_boundaries(index_col_min)), 
    //  hori_boundaries(index_row_min[index_col_min] + 1) - hori_boundaries(index_row_min[index_col_min]) + 1, 
    //  verti_boundaries(index_col_min + 1) - verti_boundaries(index_col_min) + 1));

    //-----------------------------------------  OpenCV_FAST DETECTOR ------------------------------------ //
    std::vector<cv::KeyPoint> keyPoints;
    //cv::FASTX(im_block, keyPoints, 50, true, 2);
    cv::FastFeatureDetector fast(60); // define detector threshold, TYPE_5_8 = 0, TYPE_7_12 = 1, TYPE_9_16 = 2(default)
    // feature point detection
    fast.detect(im_block, keyPoints);
    //drawKeypoints(frame.data, keyPoints, frame.data, cv::Scalar::all(255), cv::DrawMatchesFlags::DRAW_OVER_OUTIMG);
    //imshow("FAST feature", frame.data);
    //-----------------------------------------  OpenCV_FAST DETECTOR ------------------------------------ //


    //// draw keyPoints on the image
    //cv::Mat featureIm;
    //cv::drawKeypoints(im_block, keyPoints, featureIm, cv::Scalar(255, 0, 255), cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    //imshow("feature image ", featureIm);
    //cv::waitKey(0);
    //////

    num_keyPoints = keyPoints.size();
    if (num_keyPoints == 0)
      // Put a high value to use other image blocks next time(threre are no features here)
      initialized(index_col_min, index_row_min[index_col_min]) = 100;
    else
    {
      // Recover the position of the feature points in the whole image
      for (i = 0; i < num_keyPoints; ++i)
      {
        keyPoints.at(i).pt.x += (float)hori_boundaries(index_row_min[index_col_min]);
        keyPoints.at(i).pt.y += (float)verti_boundaries(index_col_min);
      }

      // For every fast corner check if there is already close features
      VectorXd no_close_features;
      no_close_features.setOnes(num_keyPoints);
      for (i = 0; i < num_keyPoints; ++i)
      {
        for (j = 0; j < uv_pred_rows; ++j)
        {
          if (sqrt(pow((double)(keyPoints.at(i).pt.x - uv_pred(j, 0)), 2) + pow((double)(keyPoints.at(i).pt.y - uv_pred(j, 1)), 2)) < 5)
            no_close_features(i) = 0;
        }
      }


      // Find how many zeros and non-zeros in the array no_close_features
      int num_close_feat = 0;
      for (i = 0; i < num_keyPoints; ++i)
      {
        if (no_close_features(i) == 0)
          num_close_feat ++;
      }

      int num_far_feat = num_keyPoints - num_close_feat;
      int *index_far_feat = new int[num_far_feat];
      int pos_close = 0, pos_far = 0;
      for (i = 0; i < num_keyPoints; ++i)
      {
        //if (no_close_features(i) == 0)
        //{
        //    index_close_feat[pos_close] = i;
        //    pos_close ++;
        //}
        if (no_close_features(i) != 0)
        {
          index_far_feat[pos_far] = i;
          pos_far ++; // not used below
        }

      }
      // If num_close_feat is equal to num_keyPoints, detection fails.
      if (num_close_feat == num_keyPoints)
        initialized(index_col_min, index_row_min[index_col_min]) = 100;

      else
      {
        // For the consideration of endoscopic image, we should eliminate the features, which are not
        // located in the circle. $$$$$$$$
        for (i = 0; i < num_far_feat; ++i)
        {
          int x = keyPoints.at(index_far_feat[i]).pt.x;
          int y = keyPoints.at(index_far_feat[i]).pt.y;
          //int cx = mono_slam->Cam->nCols_/2;
          int cx = 352;
          //int cy = mono_slam->Cam->nRows_/2;
          int cy = 255;
          //int radius = 330;
          int radius = 350;
          if (sqrt((double)((x-cx) * (x-cx) + (y-cy) * (y-cy))) > radius - 10)
          {
            for (int ii = i; ii < num_far_feat - 1; ++ii)
            {
              index_far_feat[ii]  = index_far_feat[ii + 1];
            }
            num_far_feat--; // num_far_feat is the size of the array index_far_feat
            i--; // One element has been deleted and shifted so the index should also be shifted
          }
        }
        if (num_far_feat == 0)
          initialized(index_col_min, index_row_min[index_col_min]) = 100;
        else 
        {

          //
          detected_new = true;
        // Choose a random feature
        rand_index = rand() % num_far_feat;
        index_feat_to_be_added = index_far_feat[rand_index];
        //index_feat_to_be_added = 4;
        *feature_TBadded = keyPoints.at(index_feat_to_be_added);
        // Add the feature to the filter bank
        for (i = 0; i < filter_bank_size; ++i)
        {
          add_features_inverse_depth(*feature_TBadded, &filter_bank.at(i)->x_k_k_, 
            &filter_bank.at(i)->p_k_k_, mono_slam->Cam, filter_bank.at(i)->std_z_, initial_rho, std_rho, &newFeature);
        }
        static int ind_feat = 0; // index of the features to be added
        feature_added = true;
        MonoSLAM::feature_info feat_info;
        // rowRange and colRange include right edge exclude left edge, after debugging the rowRange and colRange seem do not work properly.
        // So the subimage method is changed to cv::Rect
        //frame.data.convertTo(frame.data, CV_64FC1);
        //double pixel = frame.data.at<double>((*feature_TBadded).pt.y- half_patch_size_when_initialized, (*feature_TBadded).pt.x- half_patch_size_when_initialized);

        // Copy the subimage to the feat_info but just point to the position
        cv::Rect rect((int)((*feature_TBadded).pt.x) - half_patch_size_when_initialized, 
          (int)((*feature_TBadded).pt.y) - half_patch_size_when_initialized, 2 * half_patch_size_when_initialized + 1, 2 * half_patch_size_when_initialized + 1);
        frame.data(rect).copyTo(feat_info.patch_when_initialized);

        //feat_info.patch_when_initialized = frame.data(cv::Rect((int)((*feature_TBadded).pt.x) - half_patch_size_when_initialized, 
        //  (int)((*feature_TBadded).pt.y) - half_patch_size_when_initialized, 2 * half_patch_size_when_initialized + 1, 2 * half_patch_size_when_initialized + 1));
        //feat_info.patch_when_initialized = frame.data.
        //rowRange((int)(*feature_TBadded).pt.y - half_patch_size_when_initialized - 1, (int)(*feature_TBadded).pt.y + half_patch_size_when_initialized).
        //colRange((int)(*feature_TBadded).pt.x - half_patch_size_when_initialized - 1, (int)(*feature_TBadded).pt.x + half_patch_size_when_initialized);

        // Print out the image matrix
        //feat_info.patch_when_initialized.convertTo(feat_info.patch_when_initialized, CV_64FC1);
        //for  (int l = 0; l < 13; ++l)
        //{
        //	for (int m = 0; m < 13; ++m)
        //	{
        //		cout << (int)feat_info.patch_when_initialized.at<double>(l,m) << " ";
        //	}
        //	cout << endl;
        //}
        feat_info.ind_feat = ind_feat;
        ind_feat ++;
        feat_info.patch_when_matching = MatrixXd::Zero(2*half_patch_size_when_matching+1, 2*half_patch_size_when_matching+1);
        feat_info.r_wc_when_initialized = mono_slam->x_k_k_output.segment(5, 3);
        feat_info.R_wc_when_initialized = Q2R(mono_slam->x_k_k_output.segment(8, 4));
        feat_info.uv_when_initialized.resize(1,2);
        feat_info.uv_when_initialized << (double)(*feature_TBadded).pt.x, (double)(*feature_TBadded).pt.y;
        feat_info.half_patch_size_when_initialized = half_patch_size_when_initialized;
        feat_info.half_patch_size_when_matching = half_patch_size_when_matching;
        feat_info.times_measured = 0;
        feat_info.times_predicted = 0;
        feat_info.init_frame = step;
        feat_info.init_measurement.resize(2);
        feat_info.init_measurement << (double)(*feature_TBadded).pt.x, (double)(*feature_TBadded).pt.y;
        feat_info.type = "inversedepth";
        feat_info.yi = newFeature;
        mono_slam->features_info.push_back(feat_info);
      }
      }
    }
  }

}

void MapManagement::add_features_inverse_depth(cv::KeyPoint feature, VectorXd *x_k_k, MatrixXd *p_k_k, 
  Camera *cam, double std_z, double initial_rho, double std_rho, VectorXd *newFeature)
{
  int size_x_k_k, i;
  VectorXd Xv;
  MatrixXd P = (*p_k_k);
  size_x_k_k = (*x_k_k).size();
  Xv = (*x_k_k).head(18); // Assume that we have a 18 size of camera state
  hinv(feature, Xv, cam, initial_rho, newFeature);

  // add feature to the state vector
  (*x_k_k).conservativeResizeLike(VectorXd::Ones(size_x_k_k + 6));
  for (i = 0; i < 6; ++i)
  {
    (*x_k_k)(i + size_x_k_k) = (*newFeature)(i);
  }
  (*p_k_k) = add_a_feature_covariance_inverse_depth(P, feature, Xv, std_z, std_rho, cam);

}

void MapManagement::hinv(cv::KeyPoint feature, VectorXd Xv, Camera *cam, double initial_rho, VectorXd *newFeature)
{
  double f, U0, V0, feature_x, feature_y,
    h_LR_x, h_LR_y, h_LR_z, nx, ny, nz;
  VectorXd calibration, r_W, q_WR, h_LR(3), n(3);
  cv::KeyPoint undistort_feature;
  f = Xv(0);
  U0 = Xv(1);
  V0 = Xv(2);
  calibration = Xv.head(5);

  // Get undistorted feature points
  undistort_fm(&undistort_feature, feature, cam, calibration);
  feature_x = (double)undistort_feature.pt.x;
  feature_y = (double)undistort_feature.pt.y;

  r_W = Xv.segment(5, 3);
  q_WR = Xv.segment(8, 4);

  h_LR_x = feature_x - U0;
  h_LR_y = feature_y - V0;
  h_LR_z = f;
  h_LR << h_LR_x, h_LR_y, h_LR_z;

  n = Q2R(q_WR) * h_LR;
  nx = n(0);
  ny = n(1);
  nz = n(2);

  (*newFeature) << r_W, atan2(nx, nz), atan2(-ny, sqrt(nx * nx + nz * nz)), initial_rho; 

}

void MapManagement::undistort_fm(cv::KeyPoint *undistort_feature, cv::KeyPoint feature, Camera *cam, VectorXd calibration)
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

MatrixXd MapManagement::add_a_feature_covariance_inverse_depth(MatrixXd P, cv::KeyPoint feature, VectorXd Xv, double std_z, double std_rho, Camera *cam)
{
  double f, U0, V0, ud, vd, xd, yd, rd, uu, vu, X_c, Y_c, Z_c, 
    X_w, Y_w, Z_w;
  MatrixXd dtheta_df, dphi_df, dtheta_dCx, dphi_dCx, dtheta_dCy, dphi_dCy, 
    dtheta_dk1, dphi_dk1, dtheta_dk2, dphi_dk2;
  int cam_state_size = 18; // The size of camera state is assumed to be 18
  VectorXd calibration, q_wc, XYZ_c(3), XYZ_w(3), dtheta_dgw(3), dphi_dgw(3),
    dy_df(6), dgc_dCx(3), dy_dCx(6), dgc_dCy(3), dy_dCy(6), 
    dgc_dk1(3), dy_dk1(6), dgc_dk2(3), dy_dk2(6), dtheta_dqwr, dphi_dqwr;
  MatrixXd dgw_dqwr, dy_dqwr(6,4), dy_drw(6,3), dgw_dgc, dgc_df(3,1), dy_dxv(6,cam_state_size), 
    dyprima_dgw(5,3), dgc_dhu(3,2), dhu_dhd, dyprima_dhd, dy_dhd(6,3), R_wc,
    Ri, Padd(3,3), P_xv, P_yxv, P_y, P_xvy, P_RES, 
    P_RES00, P_RES01, P_RES02, P_RES10, P_RES11, P_RES12, P_RES20, P_RES21, P_RES22;
  cv::KeyPoint undistort_feature;

  f = Xv(0);
  U0  =  Xv(1);
  V0  =  Xv(2);
  calibration = Xv.head(5);
  q_wc = Xv.segment(8, 4);
  R_wc = Q2R(q_wc);

  ud = (double)feature.pt.x;
  vd = (double)feature.pt.y;

  xd = (ud - U0) * cam->dx_;
  yd = (vd - V0) * cam->dy_;
  rd = sqrt(xd * xd + yd * yd);

  undistort_fm(&undistort_feature, feature, cam, calibration);
  uu = (double)undistort_feature.pt.x;
  vu = (double)undistort_feature.pt.y;
  X_c = uu - U0;
  Y_c = vu - V0;
  Z_c = f;

  XYZ_c << X_c, Y_c, Z_c;

  XYZ_w = R_wc * XYZ_c;
  X_w = XYZ_w(0);
  Y_w = XYZ_w(1);
  Z_w = XYZ_w(2);

  // Derivatives
  dtheta_dgw << Z_w/(X_w * X_w + Z_w * Z_w), 0, -X_w/(X_w * X_w + Z_w * Z_w);
  dphi_dgw << (X_w * Y_w)/((X_w * X_w + Y_w * Y_w + Z_w * Z_w) * sqrt(X_w * X_w + Z_w * Z_w)),
    -sqrt(X_w * X_w + Z_w * Z_w)/(X_w * X_w + Y_w * Y_w + Z_w * Z_w), 
    (Z_w * Y_w) / ((X_w * X_w + Y_w * Y_w + Z_w * Z_w) * sqrt(X_w * X_w + Z_w * Z_w));
  dgw_dqwr = dRq_times_a_by_dq( q_wc, XYZ_c );
  //dtheta_dgw.transposeInPlace();
  //dphi_dgw.transposeInPlace();
  dtheta_dqwr = dtheta_dgw.transpose() * dgw_dqwr; // should be 1x4 vector, after arithmetic 4x1
  dphi_dqwr = dphi_dgw.transpose() * dgw_dqwr; // should be 1x4 vector, after arithmetic 4x1

  dy_dqwr << MatrixXd::Zero(3,4), dtheta_dqwr.transpose(), dphi_dqwr.transpose(), MatrixXd::Zero(1,4);
  dy_drw << MatrixXd::Identity(3,3), MatrixXd::Zero(3,3);
  dgw_dgc = R_wc;

  // Derivatives wrt f
  dgc_df << 0, 0, 1; // different from the book
  dtheta_df = dtheta_dgw.transpose() * dgw_dgc * dgc_df;
  dphi_df = dphi_dgw.transpose() * dgw_dgc * dgc_df;
  dy_df << VectorXd::Zero(3), dtheta_df, dphi_df, 0;

  // Derivatives wrt Cx
  dgc_dCx << -1, 0, 0; // different from the book
  dtheta_dCx = dtheta_dgw.transpose() * dgw_dgc * dgc_dCx;
  dphi_dCx = dphi_dgw.transpose() * dgw_dgc * dgc_dCx;
  dy_dCx << VectorXd::Zero(3), dtheta_dCx, dphi_dCx, 0;

  // Derivatives wrt Cy
  dgc_dCy << 0, -1, 0; // different from the book
  dtheta_dCy = dtheta_dgw.transpose() * dgw_dgc * dgc_dCy;
  dphi_dCy = dphi_dgw.transpose() * dgw_dgc * dgc_dCy;
  dy_dCy << VectorXd::Zero(3), dtheta_dCy, dphi_dCy, 0;

  // Derivatives wrt K1
  dgc_dk1 << (ud - U0) * pow(rd, 2), (vd - V0) * pow(rd, 2), 0; // different from the book
  dtheta_dk1 = dtheta_dgw.transpose() * dgw_dgc * dgc_dk1;
  dphi_dk1 = dphi_dgw.transpose() * dgw_dgc * dgc_dk1;
  dy_dk1 << VectorXd::Zero(3), dtheta_dk1, dphi_dk1, 0;

  // Derivatives wrt k2
  dgc_dk2 << (ud - U0) * pow(rd, 4), (vd - V0) * pow(rd, 4), 0; // different from the book
  dtheta_dk2 = dtheta_dgw.transpose() * dgw_dgc * dgc_dk2;
  dphi_dk2 = dphi_dgw.transpose() * dgw_dgc * dgc_dk2;
  dy_dk2 << VectorXd::Zero(3), dtheta_dk2, dphi_dk2, 0;

  dy_dxv << dy_df, dy_dCx, dy_dCy, dy_dk1, dy_dk2, dy_drw, dy_dqwr, MatrixXd::Zero(6,6);
  dyprima_dgw << MatrixXd::Zero(3,3), dtheta_dgw.transpose(), dphi_dgw.transpose();
  dgc_dhu << 1, 0,
    0, 1,
    0, 0;
  dhu_dhd = jacob_undistort_fm(cam, feature, calibration);
  dyprima_dhd = dyprima_dgw*dgw_dgc*dgc_dhu*dhu_dhd;
  dy_dhd << dyprima_dhd, MatrixXd::Zero(5, 1), MatrixXd::Zero(1, 2), 1;
  Ri = MatrixXd::Identity(2,2) * std_z * std_z;
  Padd << Ri, MatrixXd::Zero(2,1), MatrixXd::Zero(1,2), std_rho * std_rho;

  // Divide P into blocks
  int P_rows, P_cols;
  P_rows = P.rows(); P_cols = P.cols();
  P_xv   = P.block(0, 0, cam_state_size, cam_state_size);
  P_yxv  = P.block(cam_state_size, 0, P_rows - cam_state_size, cam_state_size);
  P_y    = P.block(cam_state_size, cam_state_size, P_rows - cam_state_size, P_cols - cam_state_size);
  P_xvy  = P.block(0, cam_state_size, cam_state_size, P_cols - cam_state_size);

  P_RES00 = P_xv; P_RES01 = P_xvy; P_RES02 = P_xv * dy_dxv.transpose();
  P_RES10 = P_yxv; P_RES11 = P_y; P_RES12 = P_yxv * dy_dxv.transpose();
  P_RES20 = dy_dxv * P_xv; P_RES21 = dy_dxv * P_xvy; 
  P_RES22 = dy_dxv * P_xv * dy_dxv.transpose() + dy_dhd * Padd * dy_dhd.transpose();

  int P_RES_size = P_RES00.cols() + P_RES01.cols() + P_RES02.cols();
  P_RES.resize(P_RES_size, P_RES_size);
  P_RES << P_RES00, P_RES01, P_RES02,
    P_RES10, P_RES11, P_RES12,
    P_RES20, P_RES21, P_RES22;
  return P_RES;

}

MatrixXd MapManagement::jacob_undistort_fm(Camera *cam, cv::KeyPoint feature, VectorXd calibration)
{
  double Cx, Cy, k1, k2, dx, dy, ud, vd, xd, yd, rd2, rd4, uu_ud, vu_vd, uu_vd, vu_ud;
  MatrixXd J_undistort(2,2);
  Cx = calibration(1);
  Cy = calibration(2);
  k1=calibration(3);
  k2=calibration(4);
  dx = cam->dx_; dy = cam->dy_;
  ud = (double)feature.pt.x;
  vd = (double)feature.pt.y;
  xd = (ud - Cx) * dx;
  yd = (vd - Cy) * dy;

  rd2 = xd * xd + yd * yd;
  rd4 = rd2 * rd2;

  uu_ud=(1+k1*rd2+k2*rd4)+(ud-Cx)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);
  vu_vd=(1+k1*rd2+k2*rd4)+(vd-Cy)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);

  uu_vd=(ud-Cx)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);
  vu_ud=(vd-Cy)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);

  J_undistort << uu_ud, uu_vd, vu_ud, vu_vd;

  return J_undistort;
}


MatrixXd MapManagement::predict_camera_measurements(VectorXd x_k_k, KalmanFilter *kalmanfilter, MonoSLAM *mono_slam)
{

  int size_x_k_k = x_k_k.size(), i, pos = 0;
  int map_features_number = mono_slam->features_info.size();
  VectorXd cali_para, t_wc, features, yi3d, yi6d, zi;
  MatrixXd r_wc;
  //mono_slam->predicted_measurements = -mono_slam->predicted_measurements.setOnes(map_features_number, 2);
  MatrixXd predicted_measurements;
  predicted_measurements = -predicted_measurements.setOnes(map_features_number, 2);

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
        predicted_measurements.block(i, 0, 1, 2) = zi.transpose();
    }
    if (mono_slam->features_info.at(i).type.compare("inversedepth") == 0)
    {
      yi6d = features.segment(pos, 6);
      pos += 6;
      if (hi_inverse_depth(&zi, yi6d, t_wc, r_wc, mono_slam->Cam, mono_slam->features_info, cali_para))
        predicted_measurements.block(i, 0, 1, 2) = zi.transpose();
    }
  }
  return predicted_measurements;
}

bool MapManagement::hi_cartesian(VectorXd *zi, VectorXd yi3d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
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

bool MapManagement::hi_inverse_depth(VectorXd *zi, VectorXd yi6d, VectorXd t_wc, MatrixXd r_wc, Camera *cam, 
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

VectorXd MapManagement::hu(VectorXd yi, Camera *cam, VectorXd cali_para)
{
  double u0, v0, f, ku, kv;
  VectorXd uv_u(2,1);
  u0 = cali_para(1);
  v0 = cali_para(2);
  f = cali_para(0);
  ku = 1/cam->dx_;
  kv = 1/cam->dy_;

  uv_u << u0 + (yi(0)/yi(2))*f, v0 + (yi(1)/yi(2))*f;
  return uv_u;
}

Vector2d MapManagement::distort_fm(VectorXd uv, Camera *cam, VectorXd cali_para)
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

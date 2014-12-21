#ifndef DATAASSOCIATOR_H
#define DATAASSOCIATOR_H
#include "Monoslam.h"

using namespace std;

class DataAssociator {
 public:
  DataAssociator();
  ~DataAssociator();

  void FilterBankMatching(MonoSLAM *mono_slam, VectorXd mu_k_km1, Frame frame);
  void Matching(Frame frame, MatrixXd predicted_measurements, MatrixXd S_predicted, 
	  vector<MonoSLAM::feature_info> features_info, Camera *cam, MatrixXd *measurements, 
	  VectorXd *measurements_map, VectorXd *predicted_measurements_map, VectorXi *matched_index);
  void JointCompatibility(MatrixXd predicted, MatrixXd measurements, MatrixXd S_predicted);


};

#endif
#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // Sum squared difference between estimation and ground truth.
  for (int i=0; i < ground_truth.size(); i++){
    // Calculate error.
    VectorXd error = ground_truth[i] - estimations[i];

    // Square the error.
    error = error.array() * error.array();

    // Sum the error.
    rmse += error;
  }

  // Find mean error.
  rmse = rmse/ground_truth.size();

  // Find square root of mean squared error.
  rmse = rmse.array().sqrt();

  return rmse;
}

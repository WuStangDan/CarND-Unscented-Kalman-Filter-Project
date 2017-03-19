#include <iostream>
#include "ukf.h"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // This will be changed after first measurement.
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Number of state variables.
  n_x_ = 5;

  // Number of augmented state variables.
  n_aug_ = 7;

  // x sigma predictions.  
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Sigma weigths.
  weights_ = VectorXd(2*n_aug_+1);
  // Spreading parameter. 
  lambda_ = 3 - n_aug_;
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if(!is_initialized_){

    previous_timestamp_ = meas_package.timestamp_;
    
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rod = meas_package.raw_measurements_[2];
      x_ << ro*cos(phi), ro*sin(phi), rod, phi, 0;
    }

    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0,0,0;
    }
    
    is_initialized_ = true;
    return;
  }

  // Compute time from last measurement in seconds.  
  float dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = meas_package.timestamp_;
  
  Prediction(dt);


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  
  - Generate sigma points,
  - Predict sigma points,
  - predict state and state covariance.
  */

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  MatrixXd x_aug = VectorXd(n_aug_);

  // Create augmented mean state.
  x_aug << x_, 0,0;

  // Create augmented covarinace matrix.
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

  // Create Square root of P matrix.
  MatrixXd Psq = P_aug.llt().matrixL();

  // Create augmented sigma points.
  Xsig_aug.col(0) = x_aug;

  for(int i=1; i<=n_aug_; i++)
  {
    Xsig_aug.col(i) = x_aug + sqrt(lambda_+n_aug_)*Psq.col(i-1);
    Xsig_aug.col(i+n_aug_) = x_aug - sqrt(lambda_+n_aug_)*Psq.col(i-1);
  }


  // Predict sigma points.
  for(int i=0; i<(2*n_aug_+1); i++)
  {
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double phi = Xsig_aug(3,i);
    double phid = Xsig_aug(4,i);
    double vak = Xsig_aug(5,i);
    double vpk = Xsig_aug(6,i);

    // Avoid division by zero.
    if (fabs(phid) > 0.001) {
      Xsig_pred_(0,i) = Xsig_aug(0,i) + v/phid*(sin(phi+phid*delta_t)-sin(phi)) + 0.5*(delta_t*delta_t)*cos(phi)*vak;
      Xsig_pred_(1,i) = Xsig_aug(1,i) + v/phid*(cos(phi)-cos(phi+phid*delta_t)) + 0.5*(delta_t*delta_t)*sin(phi)*vak;
    }
    else {
      Xsig_pred_(0,i) = Xsig_aug(0,i) + v*cos(phi)*delta_t + 0.5*(delta_t*delta_t)*cos(phi)*vak;
      Xsig_pred_(1,i) = Xsig_aug(1,i) + v*sin(phi)*delta_t + 0.5*(delta_t*delta_t)*sin(phi)*vak;
    }

    Xsig_pred_(2,i) = Xsig_aug(2,i) + delta_t*vak;
    Xsig_pred_(3,i) = Xsig_aug(3,i) + phid*delta_t+ 0.5*delta_t*delta_t*vpk;
    Xsig_pred_(4,i) = Xsig_aug(4,i) + delta_t*vpk;    
  }
  

  // Set weights.
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i=0; i<(2*n_aug_+1); i++)
  {
    weights_(i) = 1/(2*(lambda_+n_aug_));
  }

  // Predict state.
  for(int i=0; i<(2*n_aug_+1); i++)
  {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // Predict state covariance matrix.
  for(int i=0; i<(2*n_aug_+1); i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    // Normalize angles between -2pi and 2pi.
    while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
    while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  // Sigma points into measurement space.  
  for(int i=0; i<(2*n_aug_+1); i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double phi = Xsig_pred_(3,i);
    double phid = Xsig_pred_(4,i);

    Zsig(0,i) = sqrt(px*px+py*py);
    Zsig(1,i) = atan2(py,px);
    Zsig(2,i) = (px*cos(phi)*v + py*sin(phi)*v)/Zsig(0,i);  
  }

  // State measurement prediction.
  z_pred.fill(0);
  for(int i=0; i<(2*n_aug_+1); i++)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  // State covariance measurement prediction.
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;
  
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0);
  
  for(int i=0; i<(2*n_aug_+1); i++)
  {
    MatrixXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization between -2pi and 2pi.
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R;


  // Cross correlation matrix.
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i=0; i<(2*n_aug_+1); i++)
  {
    MatrixXd z_diff = Zsig.col(i) - z_pred;
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * z_diff.transpose();
  }

  // Kalman gain.
  MatrixXd K = Tc * S.inverse();


  // Measurements.
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  // Update state and state covariance.
  x_ += K*(z - z_pred);
  P_ += P_ - K*S*K.transpose();
}

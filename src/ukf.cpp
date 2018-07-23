#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);
  x_(2) = 2;

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.2;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  VectorXd v3(3);
  v3 << std_radr_ * std_radr_, std_radphi_ * std_radphi_, std_radrd_ * std_radrd_;

  R_radar_ = v3.asDiagonal();

  VectorXd v2(2);
  v2 << std_laspx_ * std_laspx_, std_laspy_ * std_laspy_;

  R_laser_ = v2.asDiagonal();
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    } else {
      float p = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);

      x_(0) = p * cos(phi);
      x_(1) = p * sin(phi);
    }

    previous_timestamp_ = meas_package.timestamp_;

    is_initialized_ = true;
    return;
  }

  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Generate augmented sigma points from state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd A = P_aug.llt().matrixL();

  // initialize sigma points
  MatrixXd Xsig = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //set first column of sigma point matrix
  Xsig.col(0) = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++) {
    Xsig.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_aug = Xsig.col(i);
    VectorXd x = x_aug.head(n_x_);
    VectorXd noise = x_aug.tail(2);

    if (fabs(x(n_x_ - 1)) > 0.00001) {
      x(0) += x(2) / x(4) * (sin(x(3) + x(4) * delta_t) - sin(x(3))) + .5 * (delta_t * delta_t) * cos(x(3)) * noise(0);
      x(1) += x(2) / x(4) * (-cos(x(3) + x(4) * delta_t) + cos(x(3))) + .5 * (delta_t * delta_t) * sin(x(3)) * noise(0);
      x(2) += delta_t * noise(0);
      x(3) += x(4) * delta_t + .5 * (delta_t * delta_t) * noise(1);
      x(4) += delta_t * noise(1);
    } else {
      x(0) += x(2) * cos(x(3)) * delta_t + .5 * (delta_t * delta_t) * cos(x(3)) * noise(0);
      x(1) += x(2) * sin(x(3)) * delta_t + .5 * (delta_t * delta_t) * sin(x(3)) * noise(0);
      x(2) += delta_t * noise(0);
      x(3) += .5 * (delta_t * delta_t) * noise(1);
      x(4) += delta_t * noise(1);
    }
    Xsig_pred_.col(i) = x;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  //calculate mean predicted measurement
  z_pred = (Zsig.array().rowwise() * weights_.transpose().array()).rowwise().sum();

  //calculate innovation covariance matrix S and cross correlation matrix Tc
  S.fill(0.0);
  S += R_laser_;

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd zdiff = Zsig.col(i) - z_pred;
    S += weights_(i) * zdiff * zdiff.transpose();

    // state difference
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (xdiff(3) > M_PI)
      xdiff(3) -= 2. * M_PI;
    while (xdiff(3) < -M_PI)
      xdiff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * xdiff * zdiff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd zdiff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * zdiff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    if (px * px + py * py < 0.0001) {
      Zsig(0, i) = 0;
      if (fabs(px) < 0.0001 && fabs(py) < 0.0001) {
        Zsig(1, i) = 0;
      } else {
        Zsig(1, i) = atan2(py, px);
      }
      Zsig(2, i) = 0;
    } else {
      Zsig(0, i) = sqrt(px * px + py * py);
      Zsig(1, i) = atan2(py, px);
      Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / Zsig(0, i);
    }
  }
  //calculate mean predicted measurement
  z_pred = (Zsig.array().rowwise() * weights_.transpose().array()).rowwise().sum();

  //calculate innovation covariance matrix S and cross correlation matrix Tc
  S.fill(0.0);
  S += R_radar_;

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd zdiff = Zsig.col(i) - z_pred;
    while (zdiff(1) > M_PI)
      zdiff(1) -= 2.0 * M_PI;
    while (zdiff(1) < -M_PI)
      zdiff(1) += 2.0 * M_PI;
    S += weights_(i) * zdiff * zdiff.transpose();

    // state difference
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (xdiff(3) > M_PI)
      xdiff(3) -= 2. * M_PI;
    while (xdiff(3) < -M_PI)
      xdiff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * xdiff * zdiff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd zdiff = z - z_pred;

  while (zdiff(1) > M_PI)
    zdiff(1) -= 2. * M_PI;
  while (zdiff(1) < -M_PI)
    zdiff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * zdiff;
  P_ = P_ - K * S * K.transpose();
}

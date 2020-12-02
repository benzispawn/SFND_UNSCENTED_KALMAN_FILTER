#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
// #include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Dimensions
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
  //       0, std_laspy_ * std_laspy_, 0, 0, 0,
  //       0, 0, 1, 0, 0,
  //       0, 0, 0, 1, 0,
  //       0, 0, 0, 0, 1;
  P_ <<     1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 0.225, 0,
            0, 0, 0, 0, 0.225;

  // Augmented dimension
  n_aug_ = n_x_ + 2;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;//30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  time_us_ = 0;

  lambda_ = 3 - n_aug_;

  size_ = 2 * n_aug_ + 1;

  weights_ = VectorXd(size_);

  // AJUST VALUE ON weights_
  double weight = 0.5/(lambda_+n_aug_);
  weights_(0) = lambda_ / (lambda_+n_aug_);

  for (int i = 1; i < size_; ++i) {
      weights_(i) = weight;
  } 

  // GENERATE SOME VETOR AND MATRIX's
  // x_aug = VectorXd(n_aug_);
  // P_aug = MatrixXd(n_aug_, n_aug_);
  // Xsig_aug = MatrixXd(n_aug_, size_);
  Xsig_pred_ = MatrixXd(n_x_, size_);
  // Q_ = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
  // Q_ << std_a_ * std_a_, 0,
  //       0, std_yawdd_ * std_yawdd_;
  
  

  radar_noise_m_ = MatrixXd(3, 3);
  radar_noise_m_ << std_radr_ * std_radr_,                         0,                       0,
                                        0, std_radphi_ * std_radphi_,                       0,
                                        0,                         0, std_radrd_ * std_radrd_;  

  lidar_noise_m_ = MatrixXd(2, 2);
  lidar_noise_m_ << std_laspx_ * std_laspx_,                       0,
                                          0, std_laspy_ * std_laspy_;  

  time_us_ = 0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
 
  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_[0], 
            meas_package.raw_measurements_[1], 
                                            0, 
                                            0, 
                                          0.0;
      // P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
      //       0, std_laspy_ * std_laspy_, 0, 0, 0,
      //       0, 0, 1, 0, 0,
      //       0, 0, 0, 1, 0,
      //       0, 0, 0, 0, 1;
      // P_ << 1, 0, 0, 0, 0,
      //       0, 1, 0, 0, 0,
      //       0, 0, 1, 0, 0,
      //       0, 0, 0, .00225, 0,
      //       0, 0, 0, 0, .00225;
                                        
    } else {
      // COORDINATE TRANSFORMATION POLAR TO RECT
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      //x_ << x, y, v, rho, rho_dot;
      x_ << x, y, v, 0, 0;

      // P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
      //   0, std_radr_ * std_radr_, 0, 0, 0,
      //   0, 0, std_radrd_ * std_radrd_, 0, 0,
      //   0, 0, 0, std_radphi_, 0,
      //   0, 0, 0, 0, std_radphi_;

      // P_ << 1, 0, 0, 0, 0,
      //       0, 1, 0, 0, 0,
      //       0, 0, 1, 0, 0,
      //       0, 0, 0, .00225, 0,
      //       0, 0, 0, 0, .00225;
    }
    //time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  //cout << meas_package.raw_measurements_ << endl;

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  // auto dt = static_cast<double>((meas_package.timestamp_ - time_us_) * 1e-6);
  time_us_ = meas_package.timestamp_;

  //  while (dt > 0.1) {
  //       constexpr double delta_t = 0.05;
  //       Prediction(delta_t);
  //       dt -= delta_t;
  //   }

  Prediction(dt);

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  // CREATE AUGMENTED MEAN VECTOR
  x_aug = VectorXd::Zero(n_aug_);

  // x_aug.fill(0.0);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // AUG STATE COVARIANCE
  P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  // CREATE SIGMA POINTS
  Xsig_aug = MatrixXd::Zero(n_aug_, size_);

  // P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // SQUARE ROOT MATRIX
  SRM_ = P_aug.llt().matrixL();

  //Xsig_aug.fill(0.0); ?????
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col( i + 1) = x_aug + sqrt(lambda_ + n_aug_) * SRM_.col(i);
    Xsig_aug.col( i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * SRM_.col(i);
  }

  // Xsig_pred_.fill(0.0);

  for (int i = 0; i < size_; i++)
  {
    // cout << Xsig_aug << endl;
    double _x = Xsig_aug(0, i);
    double _y = Xsig_aug(1, i);
    double _v = Xsig_aug(2, i);
    double _yaw = Xsig_aug(3, i);
    double _yawd = Xsig_aug(4, i);
    double _noise = Xsig_aug(5, i);
    double _noise_yaw = Xsig_aug(6, i);

    double p_x, p_y;

    if (fabs(_yawd) > 0.001)
    {
      p_x = _x + _v / _yawd * (sin(_yaw + _yawd * delta_t) - sin(_yaw));
      p_y = _y + _v / _yawd * (cos(_yaw) - cos(_yaw + _yawd * delta_t) );
    } else {
      p_x = _x + _v * delta_t * cos(_yaw);
      p_y = _y + _v * delta_t * sin(_yaw);
    }

    double v_p = _v;
    double yaw_p = _yaw + _yawd * delta_t;
    double yawd_p = _yawd;

    p_x = p_x + 0.5 * _noise * delta_t * delta_t * cos(_yaw);
    p_y = p_y +  0.5 * _noise * delta_t * delta_t * sin(_yaw);
    v_p = v_p + _noise * delta_t;
    // cout << "delta_t" << delta_t << endl;
    yaw_p = yaw_p + 0.5 * _noise_yaw * delta_t * delta_t; 
    yawd_p = yawd_p + _noise_yaw * delta_t;

    // // ADDING SOME NOISE
    Xsig_pred_(0, i) = p_x;//p_x + 0.5 * _noise * delta_t * delta_t * cos(_yaw);
    Xsig_pred_(1, i) = p_y;//p_y + 0.5 * _noise * delta_t * delta_t * sin(_yaw); 
    Xsig_pred_(2, i) = v_p;//_v + _noise * delta_t;
    Xsig_pred_(3, i) = yaw_p;//_yaw + _yawd * delta_t + 0.5 * _noise_yaw * delta_t * delta_t;
    Xsig_pred_(4, i) = yawd_p;//_yawd + _noise_yaw * delta_t;

  }

  // PREDICTION OF THE STATE MEAN
  x_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  //PREDICTION OF STATE COVARIANCE
  P_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //NORMALIZATION
    // cout << M_PI << endl;
    // cout << x_diff(3) << endl;
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  // cout << "Final of prediction" << endl;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  // cout << "Entrou na updateLidar" << endl;
  int n_z = 2;

  // H_ = Eigen::MatrixXd(n_z, n_x_);
  // H_ << 1, 0, 0, 0, 0,
  //       0, 1, 0, 0, 0; 

  // MEASUREMENT COVARIANCE MATRIX
  mea_cov_ = MatrixXd(n_z, n_z);
  mea_cov_.fill(0.0);
  // SIGMA POINT IN THE MEASUREMENT SPACE
  mea_sig_ = MatrixXd(n_z, size_);
  mea_sig_.fill(0.0);
  // PREDICTION MEASUREMENT
  mean_mea_ = VectorXd(n_z);
  mean_mea_.fill(0.0);

  mea_ = meas_package.raw_measurements_;
  for (int i = 0; i < size_; i++)
  {
    // cout << "Xsig_pred_ " 
    //           << Xsig_pred_(0, i) 
    //           << endl; 
    mea_sig_(0, i) = Xsig_pred_(0, i);
    mea_sig_(1, i) = Xsig_pred_(1, i);
  }
  //Prediction of the mean measurement
  
  for (int i = 0; i < size_; i++)
  {
    //mean_mea_ += H_.col(i) * weights_(i) * mea_sig_.col(i);
    mean_mea_ += weights_(i) * mea_sig_.col(i);
  }
  // Calculate covariance
  
  for (int i = 0; i < size_; i++)
  {
    //VectorXd mea_diff_ = H_ * mea_sig_.col(i) - mean_mea_;
    VectorXd mea_diff_ = mea_sig_.col(i) - mean_mea_;
    // ANGLE NORMALIZATION
    // while(mea_diff_(1) > M_PI) mea_diff_(1) -= 2. * M_PI;
    // while(mea_diff_(1) < -M_PI) mea_diff_(1) += 2. * M_PI; 
    mea_cov_ += weights_(i) * mea_diff_ * mea_diff_.transpose();
  }

  mea_cov_ += lidar_noise_m_; 
 
  UKF::UpdateState(mea_, n_z);

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3;

  // MEASUREMENT COVARIANCE MATRIX
  mea_cov_ = MatrixXd(n_z, n_z);

  // SIGMA POINT IN THE MEASUREMENT SPACE
  mea_sig_ = MatrixXd(n_z, size_);
  // PREDICTION MEASUREMENT
  mean_mea_ = VectorXd(n_z);

  mea_ = meas_package.raw_measurements_;


  mea_sig_.fill(0.0);
  mean_mea_.fill(0.0);
  mea_cov_.fill(0.0);

  for (int i = 0; i < size_; ++i)
  {
    double _p_x = Xsig_pred_(0, i);
    double _p_y = Xsig_pred_(1, i);
    double _p_v = Xsig_pred_(2, i);
    double _p_yaw = Xsig_pred_(3, i);

    double _p_v1 = cos(_p_yaw) * _p_v;
    double _p_v2 = sin(_p_yaw) * _p_v;

    mea_sig_(0, i) = sqrt(_p_x * _p_x + _p_y * _p_y);
    mea_sig_(1, i) = atan2(_p_y, _p_x);
    if (_p_y != 0 || _p_x != 0)
    {
      mea_sig_(2, i) = (_p_x * _p_v1 + _p_y * _p_v2) / mea_sig_(0, i);
    } else 
    {
      mea_sig_(2, i) = 0;
    }
  }

  for (int i = 0; i < size_; ++i)
  {
    mean_mea_ += mea_sig_.col(i) * weights_(i);
  }

  for (int i = 0; i < size_; ++i)
  {
    VectorXd _diff = mea_sig_.col(i) - mean_mea_;

    while (_diff(1) > M_PI) _diff(1) -= 2.*M_PI;
    while (_diff(1) < -M_PI) _diff(1) += 2.*M_PI;

    mea_cov_ += weights_(i) * _diff * _diff.transpose();
  }

  mea_cov_ += radar_noise_m_;

  UpdateState(mea_, n_z);
}

void UKF::UpdateState(const VectorXd& z, int n_z_) {
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);

  Tc_.fill(0.0);
  for (int i = 0; i < size_; i++) 
  {
    // VectorXd mean_diff = H_ * mea_sig_.col(i) - mean_mea_;
    VectorXd mean_diff = mea_sig_.col(i) - mean_mea_;

    if (n_z_ == 3) {
      while(mean_diff(1) > M_PI) mean_diff(1) -= 2.*M_PI;
      while(mean_diff(1) < -M_PI) mean_diff(1) += 2.*M_PI;
    }
    // STATE DIFFERENCE
    // VectorXd s_diff = H_ * Xsig_pred_.col(i) - x_;
    VectorXd s_diff = Xsig_pred_.col(i) - x_;

    if (n_z_ == 3) {
      while(s_diff(3) > M_PI) s_diff(3) -= 2.*M_PI;
      while(s_diff(3) < -M_PI) s_diff(3) += 2.*M_PI;
    }


    Tc_ = Tc_ + weights_(i) * s_diff * mean_diff.transpose();
  }

  // Kalman gain
  //MatrixXd T_ = MatrixXd(n_x_, n_z_);
  MatrixXd T_ = Tc_ * mea_cov_.inverse();

  //
  VectorXd r_diff = z - mean_mea_;

  //NORMALIZATION OF THE ANGLE
  if (n_z_ == 3) {
    while (r_diff(1) > M_PI) r_diff(1) -= 2.*M_PI;
    while (r_diff(1) < -M_PI) r_diff(1) += 2.*M_PI;
  }

  x_ += T_ * r_diff;
  P_ -= T_ * mea_cov_ * T_.transpose();

}
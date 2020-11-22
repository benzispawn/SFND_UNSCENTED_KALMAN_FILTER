#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
// #include <math.h>

// using Eigen::MatrixXd;
// using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // Dimensions
  n_x_ = 5;

  // initial state vector
  x_ = Eigen::VectorXd(n_x_);

  // initial covariance matrix
  P_ = Eigen::MatrixXd(n_x_, n_x_);

  // P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
  //       0, std_laspy_ * std_laspy_, 0, 0, 0,
  //       0, 0, 1, 0, 0,
  //       0, 0, 0, 1, 0,
  //       0, 0, 0, 0, 1;
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, .0225, 0,
        0, 0, 0, 0, .0225;

  // Augmented dimension
  n_aug_ = n_x_ + 2;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 10.0;//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .8;//30;
  
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

  lambda_ = 3 - n_aug_;

  size_ = 2 * n_aug_ + 1;

  weights_ = Eigen::VectorXd(size_);

  // AJUST VALUE ON weights_
  double weight = 0.5/(lambda_+n_aug_);
  weights_(0) = lambda_ / (lambda_+n_aug_);

  for (int i = 1; i < size_; ++i) {
      weights_(i) = weight;
  } 

  // GENERATE SOME VETOR AND MATRIX's
  x_aug = Eigen::VectorXd(n_aug_);
  P_aug = Eigen::MatrixXd(n_aug_, n_aug_);
  Xsig_aug = Eigen::MatrixXd(n_aug_, size_);
  Xsig_pred_ = Eigen::MatrixXd(n_x_, size_);
  Q_ = Eigen::MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
  Q_ << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;
  
  

  radar_noise_m_ = Eigen::MatrixXd(3, 3);
  radar_noise_m_ << std_radr_ * std_radr_,                         0,                       0,
                                        0, std_radphi_ * std_radphi_,                       0,
                                        0,                         0, std_radrd_ * std_radrd_;  

  lidar_noise_m_ = Eigen::MatrixXd(2, 2);
  lidar_noise_m_ << std_laspx_ * std_laspx_,                       0,
                                          0, std_laspy_ * std_laspy_;  
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
                                        
    } else {
      // COORDINATE TRANSFORMATION POLAR TO RECT
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double x = rho * std::cos(phi);
      double y = rho * std::sin(phi);
      double vx = rho_dot * std::cos(phi);
      double vy = rho_dot * std::sin(phi);
      double v = std::sqrt(vx * vx + vy * vy);
      // x_ << x, y, v, rho, rho_dot;
      x_ << x, y, v, 0, 0;

      // P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
      //   0, std_radr_ * std_radr_, 0, 0, 0,
      //   0, 0, std_radrd_ * std_radrd_, 0, 0,
      //   0, 0, 0, std_radphi_, 0,
      //   0, 0, 0, 0, std_radphi_;
    }
    is_initialized_ = true;
    return;
  }
  //std::cout << meas_package.raw_measurements_ << std::endl;

  //float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  float dt = static_cast<double>((meas_package.timestamp_ - time_us_) * 1e-6);
  time_us_ = meas_package.timestamp_;

  while (dt > 0.1)
  {
    constexpr double delta_t = 0.05;
    Prediction(dt);
    dt -= delta_t;
  }

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

void UKF::GenAugSigmaPoints(Eigen::MatrixXd& Xsig_aug)
{
/*
  // PREDICTION
  x_aug.fill(0.0);
  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // PREDICTION COVARIANCE MATRIX
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q_;
  // P_aug(n_x_, n_x_) = std_a_ * std_a_;
  // P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  // SQUARE ROOT MATRIX
  SRM_ = P_aug.llt().matrixL();

  // SIGMA POINTS
  Xsig_aug.fill(0.0);
  //std::cout << Xsig_aug << std::endl;
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + std::sqrt(lambda_ + n_aug_) * SRM_.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - std::sqrt(lambda_+n_aug_) * SRM_.col(i);
  }*/


  
}

void UKF::SigmaPointsPred(Eigen::MatrixXd& Xsig_aug, double time_t)
{
  //std::cout<< Xsig_aug < 'Depois' << std::endl;
  // PREDICTION OF SIGMA POINTS
  for (int i = 0; i < size_; i++)
  {
    double _x = Xsig_aug(0, i);
    double _y = Xsig_aug(1, i);
    double _v = Xsig_aug(2, i);
    double _yaw = Xsig_aug(3, i);
    double _yawd = Xsig_aug(4, i);
    double _noise = Xsig_aug(5, i);
    double _noise_yaw = Xsig_aug(6, i);

    double p_x, p_y;

    if (std::fabs(_yawd) > 0.001)
    {
      p_x = _x + _v / _yawd * (std::sin(_yaw + _yawd * time_t) - std::sin(_yaw));
      p_y = _y + _v / _yawd * (-1) * (std::cos(_yaw + _yawd * time_t) - std::cos(_yaw));
    } else {
      p_x = _x + _v * time_t * std::cos(_yaw);
      p_y = _y + _v * time_t * std::sin(_yaw);
    }

    double v_p = _v;
    double yaw_p = _yaw + _yawd * time_t;
    double yawd_p = _yawd;

    p_x += 0.5 * _noise * time_t * time_t * std::cos(_yaw);
    p_y += 0.5 * _noise * time_t * time_t * std::sin(_yaw);
    v_p += _noise * time_t;

    yaw_p += 0.5 * _noise_yaw * time_t * time_t; 
    yawd_p += _noise_yaw * time_t;

    // ADDING SOME NOISE
    Xsig_pred_(0, i) = p_x;//p_x + 0.5 * _noise * delta_t * delta_t * std::cos(_yaw);
    Xsig_pred_(1, i) = p_y;//p_y + 0.5 * _noise * delta_t * delta_t * std::sin(_yaw); 
    Xsig_pred_(2, i) = v_p;//_v + _noise * delta_t;
    Xsig_pred_(3, i) = yaw_p;//_yaw + _yawd * delta_t + 0.5 * _noise_yaw * delta_t * delta_t;
    Xsig_pred_(4, i) = yawd_p;//_yawd + _noise_yaw * delta_t;

  }
}

void UKF::Prediction(double ex_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  x_aug.fill(0.0);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  SRM_ = P_aug.llt().matrixL();

  int n = x_aug.size();
  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n; i++)
  {
    Xsig_aug.col( i + 1) = x_aug + std::sqrt(lambda_ + n) * SRM_.col(i);
    Xsig_aug.col( i + 1 + n) = x_aug - std::sqrt(lambda_ + n) * SRM_.col(i);
  }

  Xsig_pred_.fill(0.0);

  for (int i = 0; i < size_; i++)
  {
    double _x = Xsig_aug(0, i);
    double _y = Xsig_aug(1, i);
    double _v = Xsig_aug(2, i);
    double _yaw = Xsig_aug(3, i);
    double _yawd = Xsig_aug(4, i);
    double _noise = Xsig_aug(5, i);
    double _noise_yaw = Xsig_aug(6, i);

    double p_x, p_y;

    if (std::fabs(_yawd) > 0.001)
    {
      p_x = _x + _v / _yawd * (std::sin(_yaw + _yawd * ex_t) - std::sin(_yaw));
      p_y = _y + _v / _yawd * (-1) * (std::cos(_yaw + _yawd * ex_t) - std::cos(_yaw));
    } else {
      p_x = _x + _v * ex_t * std::cos(_yaw);
      p_y = _y + _v * ex_t * std::sin(_yaw);
    }

    double v_p = _v;
    double yaw_p = _yaw + _yawd * ex_t;
    double yawd_p = _yawd;

    p_x += 0.5 * _noise * ex_t * ex_t * std::cos(_yaw);
    p_y += 0.5 * _noise * ex_t * ex_t * std::sin(_yaw);
    v_p += _noise * ex_t;

    yaw_p += 0.5 * _noise_yaw * ex_t * ex_t; 
    yawd_p += _noise_yaw * ex_t;

    // ADDING SOME NOISE
    Xsig_pred_(0, i) = p_x;//p_x + 0.5 * _noise * delta_t * delta_t * std::cos(_yaw);
    Xsig_pred_(1, i) = p_y;//p_y + 0.5 * _noise * delta_t * delta_t * std::sin(_yaw); 
    Xsig_pred_(2, i) = v_p;//_v + _noise * delta_t;
    Xsig_pred_(3, i) = yaw_p;//_yaw + _yawd * delta_t + 0.5 * _noise_yaw * delta_t * delta_t;
    Xsig_pred_(4, i) = yawd_p;//_yawd + _noise_yaw * delta_t;

  }

  // PREDICTION OF THE STATE MEAN
  //x_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  //PREDICTION OF STATE COVARIANCE
  //P_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //NORMALIZATION
    // std::cout << M_PI << std::endl;
    //std::cout << x_diff(3) << std::endl;
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2;

  // H_ = Eigen::MatrixXd(n_z, n_x_);
  // H_ << 1, 0, 0, 0, 0,
  //       0, 1, 0, 0, 0; 

  // MEASUREMENT COVARIANCE MATRIX
  mea_cov_ = Eigen::MatrixXd(n_z, n_z);

  // SIGMA POINT IN THE MEASUREMENT SPACE
  mea_sig_ = Eigen::MatrixXd(n_z, size_);
  mea_sig_.fill(0.0);
  // PREDICTION MEASUREMENT
  mean_mea_ = Eigen::VectorXd(n_z);

  mea_ = meas_package.raw_measurements_;
  for (int i = 0; i < size_; i++)
  {
    mea_sig_(0, i) = Xsig_pred_(0, i);
    mea_sig_(1, i) = Xsig_pred_(1, i);
  }
  //Prediction of the mean measurement
  mean_mea_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    //mean_mea_ += H_.col(i) * weights_(i) * mea_sig_.col(i);
    mean_mea_ += weights_(i) * mea_sig_.col(i);
  }
  // Calculate covariance
  mea_cov_.fill(0.0);
  for (int i = 0; i < size_; i++)
  {
    //Eigen::VectorXd mea_diff_ = H_ * mea_sig_.col(i) - mean_mea_;
    Eigen::VectorXd mea_diff_ = mea_sig_.col(i) - mean_mea_;
    // ANGLE NORMALIZATION
    // while(mea_diff_(1) > M_PI) mea_diff_(1) -= 2. * M_PI;
    // while(mea_diff_(1) < -M_PI) mea_diff_(1) += 2. * M_PI; 
    mea_cov_ += weights_(i) * mea_diff_ * mea_diff_.transpose();
  }

  mea_cov_ += lidar_noise_m_; 
 
  UpdateState(mea_, 2);

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
  mea_cov_ = Eigen::MatrixXd(n_z, n_z);

  // SIGMA POINT IN THE MEASUREMENT SPACE
  mea_sig_ = Eigen::MatrixXd(n_z, size_);
  // PREDICTION MEASUREMENT
  mean_mea_ = Eigen::VectorXd(n_z);

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

    double _p_v1 = std::cos(_p_yaw) * _p_v;
    double _p_v2 = std::sin(_p_yaw) * _p_v;

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
    Eigen::VectorXd _diff = mea_sig_.col(i) - mean_mea_;

    while (_diff(1) > M_PI) _diff(1) -= 2.*M_PI;
    while (_diff(1) < -M_PI) _diff(1) += 2.*M_PI;

    mea_cov_ += weights_(i) * _diff * _diff.transpose();
  }

  mea_cov_ += radar_noise_m_;

  UpdateState(mea_, n_z);
}

void UKF::UpdateState(const Eigen::VectorXd z, int n_z_) {
  Eigen::MatrixXd Tc_ = Eigen::MatrixXd(n_x_, n_z_);

  for (int i = 0; i < size_; i++) 
  {
    // Eigen::VectorXd mean_diff = H_ * mea_sig_.col(i) - mean_mea_;
    Eigen::VectorXd mean_diff = mea_sig_.col(i) - mean_mea_;

    while(mean_diff(1) > M_PI) mean_diff(1) -= 2.*M_PI;
    while(mean_diff(1) < -M_PI) mean_diff(1) += 2.*M_PI;
    // STATE DIFFERENCE
    // Eigen::VectorXd s_diff = H_ * Xsig_pred_.col(i) - x_;
    Eigen::VectorXd s_diff = Xsig_pred_.col(i) - x_;

    while(s_diff(3) > M_PI) s_diff(3) -= 2.*M_PI;
    while(s_diff(3) < -M_PI) s_diff(3) += 2.*M_PI;

    Tc_ += weights_(i) * s_diff * mean_diff.transpose();
  }

  // Kalman gain
  Eigen::MatrixXd T_ = Eigen::MatrixXd(n_x_, n_z_);
  T_ = Tc_ * mea_cov_.inverse();

  //
  Eigen::VectorXd r_diff = z - mean_mea_;

  //NORMALIZATION OF THE ANGLE
  while (r_diff(1) > M_PI) r_diff(1) -= 2.*M_PI;
  while (r_diff(1) < -M_PI) r_diff(1) += 2.*M_PI;

  x_ += T_ * r_diff;
  P_ -= T_ * mea_cov_ * T_.transpose();

}
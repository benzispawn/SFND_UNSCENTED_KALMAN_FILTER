#ifndef UKF_H
#define UKF_H
//#define M_PI 3.141592653589799323846

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void GenAugSigmaPoints(Eigen::MatrixXd& Xsig_aug);

  void SigmaPointsPred(Eigen::MatrixXd& Xsig_aug, double delta_t);

  void UpdateState(Eigen::VectorXd z, int n_z_);

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;

  // Dimension size
  int size_;

  Eigen::VectorXd x_aug;
  Eigen::MatrixXd P_aug;
  Eigen::MatrixXd Xsig_aug;
  Eigen::MatrixXd Q_;
  //Eigen::MatrixXd H_;
  // Eigen::MatrixXd Tc_;
  // Eigen::MatrixXd T_;
  // SQUARE ROOT MATRIX
  Eigen::MatrixXd SRM_;
  //
  //Eigen::VectorXd x_diff;
  // FOR MEASUREMENTS
  Eigen::VectorXd mea_;
  Eigen::MatrixXd mea_sig_;//measurement for sigma points in measurement space
  Eigen::VectorXd mean_mea_; // mean of measurement points
  Eigen::MatrixXd mea_cov_; // covariance predicted measurement
  // NOISE MATRIX
  Eigen::MatrixXd lidar_noise_m_;
  Eigen::MatrixXd radar_noise_m_;
};

#endif  // UKF_H
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  
    x_ = F_*x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_* x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_* P_*Ht +R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_*Ht;
  MatrixXd K = PHt *Si;

  //Update state using kalman gain
  x_ = x_ + (K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1)); 
  double phi = atan2(x_(1), x_(0));
  double grad_rho = (x_(0)*x_(2) +  x_(1)*x_(3))/rho;

 
  VectorXd z_pred(3);
  z_pred << rho, phi, grad_rho;

  VectorXd y = z - z_pred;

  //Normalise the angle phi to be between -pi and +pi
  //keep subtracting pi 
  while(y(1) > M_PI || y(1) < -1.0*M_PI){
    if(y(1) > M_PI) {
        y(1) -= M_PI;
    }
    else{
        y(1) += M_PI;
    }
  }
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_* P_*Ht +R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_*Ht;
  MatrixXd K = PHt *Si;

  //Update state using kalman gain
  x_ = x_ + (K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_) * P_;


}

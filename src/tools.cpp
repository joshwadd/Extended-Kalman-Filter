#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
        std::cout << "Invalid estimate or ground truth date" << std::endl;
        return rmse;
  }

  for(unsigned int i =0; i< estimations.size(); i++){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  rmse = rmse/estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /*
   * Calculate the jacobian matrix for non-linear transistion function
   * between polarcordinates and cartesian cordinates
  
    */

    MatrixXd Hj(3,4);

    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //precompute repeated calcuations

    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    //Make sure no division by zero 
    if(fabs(c1) < 0.0001){
        cout << "Calculation of Jacobain - Error - Division by zero" << endl;
        return Hj;
    }

    //Compute the Jacobain Matrix
    Hj << (px/c2), (py/c2), 0, 0,
         -(py/c1), (px/c1), 0, 0,
         py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}

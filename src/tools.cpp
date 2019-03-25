#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    if (estimations.size()==0) {
        cout << "error: estimation vector size must be non-zero" << endl;
        return rmse;
    }
    if (estimations.size()!=ground_truth.size()) {
        cout << "error: estimation and ground truth vector sizes must be equal" << endl;
        return rmse;
    }
    for (int i=0; i < estimations.size(); ++i) {
        VectorXd diff = estimations[i] - ground_truth[i];
        rmse = rmse.array() + diff.array()*diff.array();
    }
    rmse = rmse.array()/estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
    
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    float den = px*px+py*py;
    float den_sqrt = sqrt(den);
    float den_32 = pow(den,1.5);
    if (den==0|den_sqrt==0|den_32==0)
        cout << "Error - division by zero";
    
    MatrixXd Hj(3,4);
    Hj << px/den_sqrt, py/den_sqrt, 0, 0,
    -py/den, px/den, 0, 0,
    py*(vx*py-vy*px)/den_32, px*(vy*px-vx*py)/den_32, px/den_sqrt, py/den_sqrt;
    return Hj;
}

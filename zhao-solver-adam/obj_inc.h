#ifndef OBJ_INC_H
#define OBJ_INC_H

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "util.h"

using namespace Eigen;
using namespace std;

double obj_6by6(const Vector3d& omg, const vector<vector<Event>>& events, double scale_lambda, double p_norm) {
    double fObj = 0;
    double omega_x = omg(0), omega_y = omg(1), omega_z = omg(2);
    double E = sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z);

    for (const auto& event_group : events) {
        MatrixXd A(6, event_group.size());
        
        for (size_t j = 0; j < event_group.size(); ++j) {
            double t = event_group[j].t;
            Vector3d f_j = event_group[j].bearing_vec;

            Matrix3d R_j;
            double theta = t * E;
            Matrix3d omega_cross;
            omega_cross << 0, -t*omega_z, t*omega_y,
                        t*omega_z, 0, -t*omega_x,
                        -t*omega_y, t*omega_x, 0;
            if(fabs(theta) > THRESHOLD) {
                double a = sin(theta) / theta;
                double b = (1 - cos(theta)) / (theta*theta);
                Matrix3d omega_cross_2 = omega_cross * omega_cross;
                R_j = Matrix3d::Identity() + a*omega_cross + b*omega_cross_2;
            }
            else {
                R_j = Matrix3d::Identity() + omega_cross;
            }

            Vector3d f_j_p = R_j * f_j;
            A.col(j) << t*f_j_p, f_j_p;
        }
        
        MatrixXd MM = A * A.transpose();
        SelfAdjointEigenSolver<MatrixXd> eigensolver(MM);
        double smallestEV = eigensolver.eigenvalues()(0);
        smallestEV = (smallestEV >= 0) ? smallestEV : 0;
        
        double tmp = 0;
        if (fabs(p_norm - 1.0) < THRESHOLD) {
            tmp = smallestEV;
        }
        else if (fabs(p_norm - 0.5) < THRESHOLD) {
            tmp = sqrt(smallestEV);
        }
        else {
            tmp = pow(smallestEV, p_norm);
        }
        fObj += tmp * scale_lambda;
    }
    return fObj;
}

double obj_6by6_approx(const Vector3d& omg,  int n_evt, double scale_lambda, const vector<Matrix3d>& nn, 
                                        const vector<Matrix3d>& tnn, const vector<Matrix3d>& t2nn,const vector<Matrix3d>& t3nn, 
                                        const vector<Matrix3d>& t4nn, double p_norm) {
    double fObj = 0;
    double omega_x = omg(0), omega_y = omg(1), omega_z = omg(2);

    MatrixXd M(6, 6);
    Matrix3d omega_cross;
    omega_cross<< 0, -omega_z, omega_y,
                  omega_z, 0, -omega_x,
                  -omega_y, omega_x, 0;

    for (size_t i = 0; i < n_evt; ++i) {
        M.topLeftCorner(3, 3) = t2nn[i] + omega_cross * t3nn[i] + t3nn[i] * omega_cross.transpose() + omega_cross * t4nn[i] * omega_cross.transpose();
        M.topRightCorner(3, 3) = tnn[i] + omega_cross * t2nn[i] + t2nn[i] * omega_cross.transpose() + omega_cross * t3nn[i] * omega_cross.transpose();
        M.bottomLeftCorner(3, 3) = M.topRightCorner(3, 3);
        M.bottomRightCorner(3, 3) = nn[i] + omega_cross * tnn[i] + tnn[i] * omega_cross.transpose() + omega_cross * t2nn[i] * omega_cross.transpose();

        SelfAdjointEigenSolver<MatrixXd> eigensolver(M);
        double smallestEV = eigensolver.eigenvalues()(0);
        smallestEV = (smallestEV >= 0) ? smallestEV : 0;

        double tmp = 0;
        if (fabs(p_norm - 1.0) < THRESHOLD) {
            tmp = smallestEV;
        }
        else if (fabs(p_norm - 0.5) < THRESHOLD) {
            tmp = sqrt(smallestEV);
        }
        else {
            tmp = pow(smallestEV, p_norm);
        }
        fObj += tmp * scale_lambda;
    }
    return fObj;
}

void obj_grad_inc(const Vector3d& w_init, const vector<vector<Event>>& events, 
                            double scale_lambda,  double diff_var, double p_norm, double& fObj, Vector3d& vfGrad) {
    // Compute the cost with the initial parameters
    fObj = obj_6by6(w_init, events, scale_lambda, p_norm);

    // Compute numerical gradient using central difference
    for (int i = 0; i < 3; ++i) {
        Vector3d w_plus = w_init;
        Vector3d w_minus = w_init;
        w_plus(i) += diff_var;
        w_minus(i) -= diff_var;
        
        double fObj_plus = obj_6by6(w_plus, events, scale_lambda, p_norm);
        double fObj_minus = obj_6by6(w_minus, events, scale_lambda, p_norm);
        
        vfGrad(i) = (fObj_plus - fObj_minus) / (2 * diff_var);
    }
    return;
}

void obj_grad_inc_approx(const Vector3d& w_init, const vector<vector<Event>>& events,const vector<Matrix3d>& nn,
                                        const vector<Matrix3d>& tnn, const vector<Matrix3d>& t2nn,const vector<Matrix3d>& t3nn,
                                        const vector<Matrix3d>& t4nn, double scale_lambda, double diff_var, double p_norm,
                                        double& fObj, Vector3d& vfGrad) {
    int n_evt = events.size();
    // Compute the cost with the initial parameters
    fObj = obj_6by6_approx(w_init, n_evt, scale_lambda, nn, tnn, t2nn, t3nn, t4nn, p_norm);

    // Compute numerical gradient using central difference
    for (int i = 0; i < 3; ++i) {
        Vector3d w_plus = w_init;
        Vector3d w_minus = w_init;
        w_plus(i) += diff_var;
        w_minus(i) -= diff_var;
        
        double fObj_plus = obj_6by6_approx(w_plus, n_evt, scale_lambda, nn, tnn, t2nn, t3nn, t4nn, p_norm);
        double fObj_minus = obj_6by6_approx(w_minus, n_evt, scale_lambda, nn, tnn, t2nn, t3nn, t4nn, p_norm);
        
        vfGrad(i) = (fObj_plus - fObj_minus) / (2 * diff_var);
    }
    return;
}

#endif // OBJ_INC_H
#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

#define THRESHOLD 1e-10
#define SMALL_VALUE 1e-4

struct Event {
    double t;
    Vector2d img;
    Vector3d plane_nml;
    Vector3d bearing_vec;
};

inline void expmap(double theta, double t, double omega_x, double omega_y, double omega_z, Matrix3d& R) {
    Matrix3d omega_cross;
    omega_cross << 0, -t*omega_z, t*omega_y,
                t*omega_z, 0, -t*omega_x,
                -t*omega_y, t*omega_x, 0;
    if(fabs(theta) > THRESHOLD) {
        double a = sin(theta) / theta;
        double b = (1 - cos(theta)) / (theta*theta);
        Matrix3d omega_cross_2 = omega_cross * omega_cross;
        R = Matrix3d::Identity() + a*omega_cross + b*omega_cross_2;
    }
    else {
        R = Matrix3d::Identity() + omega_cross;
    }
    return;
}

void data_compression_4th_order(const vector<vector<Event>>& events, vector<Matrix3d>& nn, vector<Matrix3d>& tnn, vector<Matrix3d>& t2nn,
                          vector<Matrix3d>& t3nn, vector<Matrix3d>& t4nn) {
    int m = events.size();
    nn.resize(m);
    tnn.resize(m);
    t2nn.resize(m);
    t3nn.resize(m);
    t4nn.resize(m);

    int k = 0;
    for (const auto& event : events) {
        Matrix3d nn_sum = Matrix3d::Zero();
        Matrix3d tnn_sum = Matrix3d::Zero();
        Matrix3d t2nn_sum = Matrix3d::Zero();
        Matrix3d t3nn_sum = Matrix3d::Zero();
        Matrix3d t4nn_sum = Matrix3d::Zero();

        for (const auto& subEvent : event) {
            double t = subEvent.t;
            Vector3d f_j = subEvent.bearing_vec;
            Matrix3d nn_matrix = f_j * f_j.transpose();
            nn_sum += nn_matrix;
            tnn_sum += t * nn_matrix;
            t2nn_sum += t * t * nn_matrix;
            t3nn_sum += t * t * t * nn_matrix;
            t4nn_sum += t * t * t * t * nn_matrix;
        }

        nn[k] = nn_sum;
        tnn[k] = tnn_sum;
        t2nn[k] = t2nn_sum;
        t3nn[k] = t3nn_sum;
        t4nn[k] = t4nn_sum;
        k++;
    }
    return;
}

void data_compression_2nd_order(const vector<vector<Event>>& events, vector<Matrix3d>& nn, vector<Matrix3d>& tnn, vector<Matrix3d>& t2nn) {
    int m = events.size();
    nn.resize(m);
    tnn.resize(m);
    t2nn.resize(m);

    int k = 0;
    for (const auto& event : events) {
        Matrix3d nn_sum = Matrix3d::Zero();
        Matrix3d tnn_sum = Matrix3d::Zero();
        Matrix3d t2nn_sum = Matrix3d::Zero();

        for (const auto& subEvent : event) {
            double t = subEvent.t;
            Vector3d plane_nml = subEvent.plane_nml;
//            Vector3d plane_nml = subEvent.plane_nml.normalized();
            Matrix3d nn_matrix = plane_nml * plane_nml.transpose();
            
            nn_sum += nn_matrix;
            tnn_sum += t * nn_matrix;
            t2nn_sum += t * t * nn_matrix;
        }

        nn[k] = nn_sum;
        tnn[k] = tnn_sum;
        t2nn[k] = t2nn_sum;
        k++;
    }
    return;
}

#endif // UTIL_H
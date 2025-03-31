#ifndef OPTIMIZER_TRANSLATION_H
#define OPTIMIZER_TRANSLATION_H

#include <limits>
#include <numeric>
#include <iostream>
#include "util.h"

using namespace Eigen;
using namespace std;

bool translation_cop_core(vector<vector<Event>>& events, vector<vector<Matrix3d>>& orientations, Vector3d& v_sol, vector<Matrix3d>& line_struct_all) {
    int n = events.size();
    MatrixXd CC(n, 3);
    line_struct_all.resize(n);
    for (size_t i = 0; i < n; ++i) {
        vector<Event> event_group = events[i];
        int m = event_group.size();
        MatrixXd AA(m, 3);
        MatrixXd BB(m, 6);
        // line direction d
        for (size_t j = 0; j < m; ++j) {
            double t = event_group[j].t;
            Vector3d plane_nml = event_group[j].plane_nml;
//            Vector3d plane_nml = event_group[j].plane_nml.normalized();

            Matrix3d R_j = orientations[i][j];
            Vector3d n_j_p = R_j * plane_nml;

            AA.row(j) << n_j_p.transpose();
        }
        Eigen::JacobiSVD<MatrixXd> svdA(AA, ComputeFullV);
        Vector3d d = svdA.matrixV().col(2).normalized();

        // translation v and line moment mmt
        for (size_t j = 0; j < m; ++j) {
            double t = event_group[j].t;
            Vector3d f_j = event_group[j].bearing_vec;
            Matrix3d R_j = orientations[i][j];
            Vector3d f_j_p = R_j * f_j;
            BB.row(j) << t*f_j_p.cross(d).transpose(), f_j_p.transpose();
        }
        Eigen::JacobiSVD<MatrixXd> svdB(BB, ComputeFullV);

        Vector3d v, mmt;
        double s = svdB.matrixV().col(5).tail(3).norm();
        if ( s > SMALL_VALUE ) {
            v = svdB.matrixV().col(5).head(3);
            mmt = svdB.matrixV().col(5).tail(3);
        }
        else {
            v = svdB.matrixV().col(4).head(3) + svdB.matrixV().col(5).head(3);
            mmt = svdB.matrixV().col(4).tail(3) + svdB.matrixV().col(5).tail(3);
        }
        s = mmt.norm();
        v = v / s;
        mmt = mmt / s;
        // line-dependent reference R_l
        Vector3d e1_l = d;
        Vector3d e2_l = -mmt;
        Vector3d e3_l = e1_l.cross(e2_l);
        Matrix3d R_l;
        R_l.col(0) = e1_l;
        R_l.col(1) = e2_l;
        R_l.col(2) = e3_l;
        Vector3d u_l = R_l.transpose() * v;
        Vector3d tmp = u_l(1) * e3_l - u_l(2) * e2_l;
        CC.row(i) = tmp.transpose();

        line_struct_all[i] = R_l;
    }
    Eigen::JacobiSVD<MatrixXd> svdC(CC, ComputeFullV);
    v_sol = svdC.matrixV().col(2).normalized();

    bool is_success = true;
    if (std::isnan(v_sol(0)) || std::isnan(v_sol(1)) || std::isnan(v_sol(2))) {
        is_success = false;
    }
    return is_success;
}

bool compute_translation_cop(vector<vector<Event>>& events, Vector3d& omg_sol, Vector3d& v_sol, vector<Matrix3d>& line_struct_all) {
    int n = events.size();
    if (n < 2) {
        std::cout << "error: at least two lines are needed!" << std::endl;
        return false; 
    }
    if (std::isnan(omg_sol(0)) || std::isnan(omg_sol(1)) || std::isnan(omg_sol(2))) {
        std::cout << "error: angular velocity is NAN!" << std::endl;
        return false; 
    }

    double omega_x = omg_sol(0), omega_y = omg_sol(1), omega_z = omg_sol(2);
    double E = sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z);
    vector<vector<Matrix3d>> orientations(n);
    for (size_t i = 0; i < n; ++i) {
        vector<Event> event_group = events[i];
        int m = event_group.size();
        vector<Matrix3d> R_vec(m);
        for (size_t j = 0; j < m; ++j) {
            double t = event_group[j].t;
            Matrix3d R;
            double theta = t * E;
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
            R_vec[j] = R;
        }
        orientations[i] = R_vec;
    }

    bool is_success = translation_cop_core(events, orientations, v_sol, line_struct_all);

    return is_success;
}

bool translation_cop_core(vector<vector<Event>>& events, vector<vector<Matrix3d>>& orientations, Vector3d& v_sol)
{
    vector<Matrix3d> line_struct_all;
    return translation_cop_core(events, orientations, v_sol, line_struct_all);
}

bool compute_translation_cop(vector<vector<Event>>& events, Vector3d& omg_sol, Vector3d& v_sol)
{
    vector<Matrix3d> line_struct_all;
    return compute_translation_cop(events, omg_sol, v_sol, line_struct_all);
}

#endif // OPTIMIZER_TRANSLATION_H
/*Some of the code is derived from OpenGV, authored by Laurent Kneip.
For more details, please refer to https://github.com/laurentkneip/opengv/blob/master/src/relative_pose/modules/eigensolver/modules.cpp
*/
#ifndef OBJ_COP_FLOW_H
#define OBJ_COP_FLOW_H

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "util.h"

using namespace Eigen;
using namespace std;

void smallest_eigenvalue_grad(const Matrix3d& M, const Matrix3d& M_jac1, const Matrix3d& M_jac2, const Matrix3d& M_jac3,
                                 double& smallestEV, double& smallestEV_jac1, double& smallestEV_jac2, double& smallestEV_jac3) {
    double b = -M(0,0)-M(1,1)-M(2,2);
    double b_jac1 = -M_jac1(0,0)-M_jac1(1,1)-M_jac1(2,2);
    double b_jac2 = -M_jac2(0,0)-M_jac2(1,1)-M_jac2(2,2);
    double b_jac3 = -M_jac3(0,0)-M_jac3(1,1)-M_jac3(2,2);
    double c = -pow(M(0,2),2)-pow(M(1,2),2)-pow(M(0,1),2)+M(0,0)*M(1,1)+M(0,0)*M(2,2)+M(1,1)*M(2,2);
    double c_jac1 = -2.0*M(0,2)*M_jac1(0,2)-2.0*M(1,2)*M_jac1(1,2)-2.0*M(0,1)*M_jac1(0,1)
        +M_jac1(0,0)*M(1,1)+M(0,0)*M_jac1(1,1)+M_jac1(0,0)*M(2,2)
        +M(0,0)*M_jac1(2,2)+M_jac1(1,1)*M(2,2)+M(1,1)*M_jac1(2,2);
    double c_jac2 = -2.0*M(0,2)*M_jac2(0,2)-2.0*M(1,2)*M_jac2(1,2)-2.0*M(0,1)*M_jac2(0,1)
        +M_jac2(0,0)*M(1,1)+M(0,0)*M_jac2(1,1)+M_jac2(0,0)*M(2,2)
        +M(0,0)*M_jac2(2,2)+M_jac2(1,1)*M(2,2)+M(1,1)*M_jac2(2,2);
    double c_jac3 = -2.0*M(0,2)*M_jac3(0,2)-2.0*M(1,2)*M_jac3(1,2)-2.0*M(0,1)*M_jac3(0,1)
        +M_jac3(0,0)*M(1,1)+M(0,0)*M_jac3(1,1)+M_jac3(0,0)*M(2,2)
        +M(0,0)*M_jac3(2,2)+M_jac3(1,1)*M(2,2)+M(1,1)*M_jac3(2,2);
    double d = M(1,1)*pow(M(0,2),2)+M(0,0)*pow(M(1,2),2)+M(2,2)*pow(M(0,1),2)-
        M(0,0)*M(1,1)*M(2,2)-2*M(0,1)*M(1,2)*M(0,2);
    double d_jac1 = M_jac1(1,1)*pow(M(0,2),2)+M(1,1)*2*M(0,2)*M_jac1(0,2)
        +M_jac1(0,0)*pow(M(1,2),2)+M(0,0)*2.0*M(1,2)*M_jac1(1,2)
        +M_jac1(2,2)*pow(M(0,1),2)+M(2,2)*2.0*M(0,1)*M_jac1(0,1)
        -M_jac1(0,0)*M(1,1)*M(2,2)-M(0,0)*M_jac1(1,1)*M(2,2)
        -M(0,0)*M(1,1)*M_jac1(2,2)-2.0*(M_jac1(0,1)*M(1,2)*M(0,2)
        +M(0,1)*M_jac1(1,2)*M(0,2)+M(0,1)*M(1,2)*M_jac1(0,2));
    double d_jac2 = M_jac2(1,1)*pow(M(0,2),2)+M(1,1)*2*M(0,2)*M_jac2(0,2)
        +M_jac2(0,0)*pow(M(1,2),2)+M(0,0)*2.0*M(1,2)*M_jac2(1,2)
        +M_jac2(2,2)*pow(M(0,1),2)+M(2,2)*2.0*M(0,1)*M_jac2(0,1)
        -M_jac2(0,0)*M(1,1)*M(2,2)-M(0,0)*M_jac2(1,1)*M(2,2)
        -M(0,0)*M(1,1)*M_jac2(2,2)-2.0*(M_jac2(0,1)*M(1,2)*M(0,2)
        +M(0,1)*M_jac2(1,2)*M(0,2)+M(0,1)*M(1,2)*M_jac2(0,2));
    double d_jac3 = M_jac3(1,1)*pow(M(0,2),2)+M(1,1)*2*M(0,2)*M_jac3(0,2)
        +M_jac3(0,0)*pow(M(1,2),2)+M(0,0)*2.0*M(1,2)*M_jac3(1,2)
        +M_jac3(2,2)*pow(M(0,1),2)+M(2,2)*2.0*M(0,1)*M_jac3(0,1)
        -M_jac3(0,0)*M(1,1)*M(2,2)-M(0,0)*M_jac3(1,1)*M(2,2)
        -M(0,0)*M(1,1)*M_jac3(2,2)-2.0*(M_jac3(0,1)*M(1,2)*M(0,2)
        +M(0,1)*M_jac3(1,2)*M(0,2)+M(0,1)*M(1,2)*M_jac3(0,2));
    double s = 2*pow(b,3)-9*b*c+27*d;
    double t = 4*pow((pow(b,2)-3*c),3);
    double s_jac1 = 2.0*3.0*pow(b,2)*b_jac1-9.0*b_jac1*c-9.0*b*c_jac1+27.0*d_jac1;
    double s_jac2 = 2.0*3.0*pow(b,2)*b_jac2-9.0*b_jac2*c-9.0*b*c_jac2+27.0*d_jac2;
    double s_jac3 = 2.0*3.0*pow(b,2)*b_jac3-9.0*b_jac3*c-9.0*b*c_jac3+27.0*d_jac3;
    double t_jac1 = 4.0*3.0*pow((pow(b,2)-3.0*c),2)*(2.0*b*b_jac1-3.0*c_jac1);
    double t_jac2 = 4.0*3.0*pow((pow(b,2)-3.0*c),2)*(2.0*b*b_jac2-3.0*c_jac2);
    double t_jac3 = 4.0*3.0*pow((pow(b,2)-3.0*c),2)*(2.0*b*b_jac3-3.0*c_jac3);

    double alpha = acos(s/sqrt(t));
    double alpha_jac1 = -1.0/sqrt(1.0-(pow(s,2)/t)) * (s_jac1*sqrt(t)-s*0.5*pow(t,-0.5)*t_jac1)/t;
    double alpha_jac2 = -1.0/sqrt(1.0-(pow(s,2)/t)) * (s_jac2*sqrt(t)-s*0.5*pow(t,-0.5)*t_jac2)/t;
    double alpha_jac3 = -1.0/sqrt(1.0-(pow(s,2)/t)) * (s_jac3*sqrt(t)-s*0.5*pow(t,-0.5)*t_jac3)/t;
    double beta = alpha/3;
    double beta_jac1 = alpha_jac1/3.0;
    double beta_jac2 = alpha_jac2/3.0;
    double beta_jac3 = alpha_jac3/3.0;
    double y = cos(beta);
    double y_jac1 = -sin(beta)*beta_jac1;
    double y_jac2 = -sin(beta)*beta_jac2;
    double y_jac3 = -sin(beta)*beta_jac3;

    double r = 0.5*sqrt(t);
    double r_jac1 = 0.25*pow(t,-0.5)*t_jac1;
    double r_jac2 = 0.25*pow(t,-0.5)*t_jac2;
    double r_jac3 = 0.25*pow(t,-0.5)*t_jac3;
    double w = pow(r,(1.0/3.0));
    double w_jac1 = (1.0/3.0)*pow(r,-2.0/3.0)*r_jac1;
    double w_jac2 = (1.0/3.0)*pow(r,-2.0/3.0)*r_jac2;
    double w_jac3 = (1.0/3.0)*pow(r,-2.0/3.0)*r_jac3;

    double k = w*y;
    double k_jac1 = w_jac1*y+w*y_jac1;
    double k_jac2 = w_jac2*y+w*y_jac2;
    double k_jac3 = w_jac3*y+w*y_jac3;
    smallestEV = (-b-2*k)/3;
    smallestEV_jac1 = (-b_jac1-2.0*k_jac1)/3.0;
    smallestEV_jac2 = (-b_jac2-2.0*k_jac2)/3.0;
    smallestEV_jac3 = (-b_jac3-2.0*k_jac3)/3.0;

    smallestEV = (smallestEV >= 0) ? smallestEV : 0;

    return;
}

void obj_grad_cop_flow(const Vector3d& w_init, const vector<vector<Event>>& events, double scale_lambda, double& fObj, Vector3d& vfGrad) {
    double omega_x = w_init(0);
    double omega_y = w_init(1);
    double omega_z = w_init(2);

    double E = sqrt(omega_x * omega_x + omega_y * omega_y + omega_z * omega_z);

    fObj = 0;
    vfGrad = Vector3d::Zero();

    for (const auto& event : events) {
        Matrix3d M = Matrix3d::Zero();
        Matrix3d M_jac1 = Matrix3d::Zero(), M_jac2 = Matrix3d::Zero(), M_jac3 = Matrix3d::Zero();

        for (const auto& subEvent : event) {
            double t = subEvent.t;
            Vector3d plane_nml = subEvent.plane_nml;
//            Vector3d plane_nml = subEvent.plane_nml.normalized();

            Matrix3d nn = plane_nml * plane_nml.transpose();

            double theta = t * E;
            Matrix3d omega_cross;
            omega_cross << 0, -omega_z * t, omega_y * t,
                           omega_z * t, 0, -omega_x * t,
                           -omega_y * t, omega_x * t, 0;

            double sinTheta = sin(theta);
            double cosTheta = cos(theta);
            Matrix3d R;
            if(fabs(theta) > THRESHOLD) {
                R = Matrix3d::Identity() + (sinTheta / theta) * omega_cross + ((1 - cosTheta) / (theta * theta)) * (omega_cross * omega_cross);
            }
            else {
                R = Matrix3d::Identity() + omega_cross;
            }
            M += R * nn * R.transpose();

            // Compute derivatives
            double E2 = E * E;
            double E3 = E2 * E;
            double tE = t * E;
            double t2 = t * t;
            double cos_tE = cos(tE);
            double sin_tE = sin(tE);
            double common_denominator = t * E3;
            double dA_common_term = cos_tE / E2 - sin_tE / common_denominator;
            double dA_dx = omega_x * dA_common_term;
            double dA_dy = omega_y * dA_common_term;
            double dA_dz = omega_z * dA_common_term;

            Matrix3d dB_dx, dB_dy, dB_dz;
            dB_dx << 0, 0, 0, 0, 0, -t, 0, t, 0;
            dB_dy << 0, 0, t, 0, 0, 0, -t, 0, 0;
            dB_dz << 0, -t, 0, t, 0, 0, 0, 0, 0;
            double common_term1 = -2 * (1 - cos_tE) / (t2 * E3);
            double common_term2 = sin_tE / common_denominator;
            double dC_dx = omega_x * (common_term1 + common_term2);
            double dC_dy = omega_y * (common_term1 + common_term2);
            double dC_dz = omega_z * (common_term1 + common_term2);

            Matrix3d dB2_dx, dB2_dy, dB2_dz;
            dB2_dx << 0, t * t * omega_y, t * t * omega_z, t * t * omega_y, -2 * t * t * omega_x, 0, t * t * omega_z, 0, -2 * t * t * omega_x;
            dB2_dy << -2 * t * t * omega_y, t * t * omega_x, 0, t * t * omega_x, 0, t * t * omega_z, 0, t * t * omega_z, -2 * t * t * omega_y;
            dB2_dz << -2 * t * t * omega_z, 0, t * t * omega_x, 0, -2 * t * t * omega_z, t * t * omega_y, t * t * omega_x, t * t * omega_y, 0;
            Matrix3d omega_cross_2 = omega_cross * omega_cross;

            Matrix3d R_transpose = R.transpose();
            double sinTheta_over_theta = sinTheta / theta;
            double cosTheta_over_theta_squared = (1 - cosTheta) / (theta * theta);

            Matrix3d dR_dx = dA_dx * omega_cross + (sinTheta_over_theta) * dB_dx + dC_dx * (omega_cross_2) + (cosTheta_over_theta_squared) * dB2_dx;
            Matrix3d dR_dy = dA_dy * omega_cross + (sinTheta_over_theta) * dB_dy + dC_dy * (omega_cross_2) + (cosTheta_over_theta_squared) * dB2_dy;
            Matrix3d dR_dz = dA_dz * omega_cross + (sinTheta_over_theta) * dB_dz + dC_dz * (omega_cross_2) + (cosTheta_over_theta_squared) * dB2_dz;

            M_jac1 += dR_dx * nn * R_transpose + R * nn * dR_dx.transpose();
            M_jac2 += dR_dy * nn * R_transpose + R * nn * dR_dy.transpose();
            M_jac3 += dR_dz * nn * R_transpose + R * nn * dR_dz.transpose();
        }

        double smallestEV, smallestEV_jac1, smallestEV_jac2, smallestEV_jac3;
        smallest_eigenvalue_grad(M, M_jac1, M_jac2, M_jac3, smallestEV, smallestEV_jac1, smallestEV_jac2, smallestEV_jac3);

        fObj += smallestEV * scale_lambda;
        vfGrad(0) += smallestEV_jac1 * scale_lambda;
        vfGrad(1) += smallestEV_jac2 * scale_lambda;
        vfGrad(2) += smallestEV_jac3 * scale_lambda;
    }
    return;
}

void obj_grad_cop_flow_approx(const Vector3d& w_init, const vector<vector<Event>>& events,
                                    const vector<Matrix3d>& nn, const vector<Matrix3d>& tnn, const vector<Matrix3d>& t2nn,
                                    double scale_lambda, double& fObj, Vector3d& vfGrad) {
    double omega_x = w_init(0);
    double omega_y = w_init(1);
    double omega_z = w_init(2);

    Matrix3d omega_cross, omega_cross_dx, omega_cross_dy, omega_cross_dz;
    omega_cross << 0, -omega_z, omega_y, omega_z, 0, -omega_x, -omega_y, omega_x, 0;
    omega_cross_dx << 0, 0, 0, 0, 0, -1, 0, 1, 0;
    omega_cross_dy << 0, 0, 1, 0, 0, 0, -1, 0, 0;
    omega_cross_dz << 0, -1, 0, 1, 0, 0, 0, 0, 0;

    fObj = 0;
    vfGrad = Vector3d::Zero();

    for (size_t i = 0; i < events.size(); ++i) {
        Matrix3d M = nn[i] + omega_cross * tnn[i] + tnn[i] * omega_cross.transpose() + omega_cross * t2nn[i] * omega_cross.transpose();
        Matrix3d M_jac1 = omega_cross_dx * tnn[i] + tnn[i] * omega_cross_dx.transpose() + omega_cross_dx * t2nn[i] * omega_cross.transpose() + omega_cross * t2nn[i] * omega_cross_dx.transpose();
        Matrix3d M_jac2 = omega_cross_dy * tnn[i] + tnn[i] * omega_cross_dy.transpose() + omega_cross_dy * t2nn[i] * omega_cross.transpose() + omega_cross * t2nn[i] * omega_cross_dy.transpose();
        Matrix3d M_jac3 = omega_cross_dz * tnn[i] + tnn[i] * omega_cross_dz.transpose() + omega_cross_dz * t2nn[i] * omega_cross.transpose() + omega_cross * t2nn[i] * omega_cross_dz.transpose();

        double smallestEV, smallestEV_jac1, smallestEV_jac2, smallestEV_jac3;
        smallest_eigenvalue_grad(M, M_jac1, M_jac2, M_jac3, smallestEV, smallestEV_jac1, smallestEV_jac2, smallestEV_jac3);

        fObj += smallestEV * scale_lambda;
        vfGrad(0) += smallestEV_jac1 * scale_lambda;
        vfGrad(1) += smallestEV_jac2 * scale_lambda;
        vfGrad(2) += smallestEV_jac3 * scale_lambda;
    }
    return;
}

#endif // OBJ_COP_FLOW_H
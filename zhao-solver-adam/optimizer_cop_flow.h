#ifndef OPTIMIZER_COP_FLOW_H
#define OPTIMIZER_COP_FLOW_H

#include <limits>
#include <numeric>
#include "obj_cop_flow.h"

using namespace Eigen;
using namespace std;

double adam_cop_flow(vector<vector<Event>>& events,
        Vector3d& x, int method_type, double scale_lambda, double diff_var, double p_norm, 
        double TolFun, double TolX, int MaxIter, int MaxFunEvals, double stepSize,
        double beta1, double beta2, double epsilon, int nEpochSize)
{
    int numberOfVariables = 3;
    // Add perturbation if using exact rotation and initial values are close to zero
    if (method_type == 301) {
        if (x.array().abs().maxCoeff() < THRESHOLD) {
            x += Vector3d::Random() * SMALL_VALUE;
        }
    }

    VectorXd vfCost(MaxIter + 1);
    MatrixXd xHist(numberOfVariables, MaxIter + 1);
    xHist.col(0) = x;

    // Initialize variables for optimization
    double fval;
    Vector3d gradient;
    int funccount = 2;

    // Compute compressed data if using first-order approximation
    vector<Matrix3d> nn, tnn, t2nn;
    if (method_type == 302 || method_type == 303) {
        data_compression_2nd_order(events, nn, tnn, t2nn);
    }

    // Function to perform one iteration of optimization
    auto performIteration = [&](const auto& derivation_func) {
        // Update biased first moment estimate
        Vector3d m = Vector3d::Zero();
        // Update biased second raw moment estimate
        Vector3d v = Vector3d::Zero();

        for (int nIter = 0; nIter < MaxIter; ++nIter) {
            // Update moment estimates
            m = beta1 * m + (1 - beta1) * gradient;
            v = beta2 * v + (1 - beta2) * gradient.array().square().matrix();

            // Compute bias-corrected estimates
            double tmp = (1 - pow(beta1, nIter + 1));
            Vector3d mHat = m / tmp;
            Vector3d vHat = v / tmp;

            // Determine step to take at this iteration
            Vector3d vfStep = stepSize * mHat.array() / (vHat.array().sqrt() + epsilon);
            x -= vfStep;

            // Call derivation function
            derivation_func(x, events, scale_lambda, fval, gradient);

            vfCost(nIter + 1) = fval;
            funccount += 2;
            xHist.col(nIter + 1) = x;

            int nFirstCost = max(1, nIter + 1 - nEpochSize);
            double fEstCost = vfCost.segment(nFirstCost - 1, nIter + 2 - nFirstCost).mean();
            double fImprEst = std::abs(fEstCost - vfCost(nFirstCost - 1));

            // Check for stopping conditions
            if (nIter >= nEpochSize) {
                if (fImprEst < TolFun / nEpochSize || vfStep.norm() < TolX || funccount > MaxFunEvals) {
                    break;
                }
            }
        }
    };

    // Perform optimization based on method type
    if (method_type == 302 || method_type == 303) {
        obj_grad_cop_flow_approx(x, events, nn, tnn, t2nn, scale_lambda, fval, gradient);
        vfCost(0) = fval;
        performIteration([&](const Vector3d& x, const vector<vector<Event>>& evts, double scale, double& fval, Vector3d& grad) {
            obj_grad_cop_flow_approx(x, evts, nn, tnn, t2nn, scale, fval, grad);
        });
    }
    if (method_type == 301 || method_type == 303) {
        obj_grad_cop_flow(x, events, scale_lambda, fval, gradient);
        vfCost(0) = fval;
        performIteration([&](const Vector3d& x, const vector<vector<Event>>& evts, double scale, double& fval, Vector3d& grad) {
            obj_grad_cop_flow(x, evts, scale, fval, grad);
        });
    }
//    Vector3d vfGrad_tmp;
//    obj_grad_cop_flow(x, events, scale_lambda, objValue, vfGrad_tmp);
//    objValue /= scale_lambda;
    double objValue = fval / scale_lambda;
	
    return objValue;
}

#endif // OPTIMIZER_COP_FLOW_H
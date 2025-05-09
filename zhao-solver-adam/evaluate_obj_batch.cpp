#include <iostream>
#include "mex.h"
#include "matrix.h"
#include "obj_cop_flow_fdm.h"
#include "obj_cop_flow.h"
#include "obj_inc.h"

using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // 1st input parameter
    const mxArray *matlab_events;
    if (nrhs > 0 && !mxIsEmpty(prhs[0])) {
        matlab_events = prhs[0];
    }
    else {
        mexErrMsgIdAndTxt("MATLAB:evaluate_obj:invalidFirstInput", "The first input cannot be empty.");
    }
    
    // 2nd input parameter
    // intial value of wx for the angular velocity
    mwSize nx = mxGetNumberOfElements(prhs[1]);
    double *inputVectorX = mxGetPr(prhs[1]);
    double* wx = new double[nx];
    for (mwSize i = 0; i < nx; ++i) {
        wx[i] = inputVectorX[i];
    }

    // 3rd input parameter
    // intial value of wy for the angular velocity
    mwSize ny = mxGetNumberOfElements(prhs[2]);
    double *inputVectorY = mxGetPr(prhs[2]);
    double* wy = new double[ny];
    for (mwSize i = 0; i < ny; ++i) {
        wy[i] = inputVectorY[i];
    }

    // 4th input parameter
    // intial value of wz for the angular velocity
    mwSize nz = mxGetNumberOfElements(prhs[3]);
    double *inputVectorZ = mxGetPr(prhs[3]);
    double* wz = new double[nz];
    for (mwSize i = 0; i < nz; ++i) {
        wz[i] = inputVectorZ[i];
    }

    // 5th input parameter
    // method type
    int method_type = 103;
    if (nrhs > 4 && !mxIsEmpty(prhs[4])) {
        method_type = static_cast<int>(mxGetScalar(prhs[4]));
    }

    // 6th input parameter
    // hyperparameters for problems and optimiziers
    double scale_lambda, diff_var, p_norm;
    scale_lambda = 1e6;
    diff_var = 1e-6; // Small deviation for numerical gradient computation
    p_norm = 1.0;
    if (nrhs > 5 && !mxIsEmpty(prhs[5])) {
        const mxArray* paramsStruct = prhs[5];
        if (mxGetField(paramsStruct, 0, "scale_lambda") != NULL) {
            scale_lambda = mxGetScalar(mxGetField(paramsStruct, 0, "scale_lambda"));
        }
        if (mxGetField(paramsStruct, 0, "diff_var") != NULL) {
            diff_var = mxGetScalar(mxGetField(paramsStruct, 0, "diff_var"));
        }
        if (mxGetField(paramsStruct, 0, "p_norm") != NULL) {
            p_norm = mxGetScalar(mxGetField(paramsStruct, 0, "p_norm"));
            if (abs(p_norm-1.0) > THRESHOLD && abs(p_norm-2.0) > THRESHOLD && abs(p_norm - 0.5) > THRESHOLD) {
                mexErrMsgIdAndTxt("MATLAB:relpose_event:invalidNorm", "The supported norms are 1, 2, and 1/2.");
            }
        }
    }

    // format conversion for events
    vector<vector<Event>> cpp_events;
    for (int i = 0; i < mxGetNumberOfElements(matlab_events); ++i) {
        mxArray *mxSubEvent = mxGetCell(matlab_events, i);
        vector<Event> subEvents;
        for (int j = 0; j < mxGetNumberOfElements(mxSubEvent); ++j) {
            mxArray *subEvent = mxGetCell(mxSubEvent, j);
            Event e;
            e.t = mxGetScalar(mxGetField(subEvent, 0, "t"));
            double *img_ptr = mxGetPr(mxGetField(subEvent, 0, "img"));
            Vector2d tmp1(img_ptr[0], img_ptr[1]);
            e.img << tmp1;
            Vector3d tmp2;
            tmp2 << tmp1, 1.0;
            e.bearing_vec = tmp2.normalized();
            double *plane_nml_ptr = mxGetPr(mxGetField(subEvent, 0, "plane_nml"));
            Vector3d tmp3(plane_nml_ptr[0], plane_nml_ptr[1], plane_nml_ptr[2]);
            e.plane_nml << tmp3.normalized();
            subEvents.push_back(e);
        }
        cpp_events.push_back(subEvents);
    }

    // Output the objective
    if (nlhs == 0) {
        return;
    }

    mwSize dims[3] = {ny, nx, nz};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *arrayData = mxGetPr(plhs[0]);
    Vector3d x = Vector3d::Zero();
    double objValue = -1.0;
    Vector3d vfGrad_tmp;
    if (method_type >= 101 && method_type <= 103) {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                for (int k = 0; k < nz; k++) {
                    x << wx[j], wy[i], wz[k];
                    objValue = obj_6by6(x, cpp_events, scale_lambda, p_norm);
                    arrayData[k*nx*ny + j*ny + i] = objValue / scale_lambda;;
                }
            }
        }
    }

    if (method_type >= 301 && method_type <= 303) {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                for (int k = 0; k < nz; k++) {
                    x << wx[j], wy[i], wz[k];
                    obj_grad_cop_flow(x, cpp_events, scale_lambda, objValue, vfGrad_tmp);
                    arrayData[k*nx*ny + j*ny + i] = objValue / scale_lambda;;
                }
            }
        }
    }

    if (method_type >= 304 && method_type <= 306) {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                for (int k = 0; k < nz; k++) {
                    x << wx[j], wy[i], wz[k];
                    objValue = obj_3by3_fdm(x, cpp_events, scale_lambda, p_norm);
                    arrayData[k*nx*ny + j*ny + i] = objValue / scale_lambda;;
                }
            }
        }
    }

    return;
}
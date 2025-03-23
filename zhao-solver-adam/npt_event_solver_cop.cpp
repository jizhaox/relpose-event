#include <iostream>
#include <chrono>
#include "mex.h"
#include "matrix.h"
#include "optimizer_translation.h"

using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 0) {
        mexPrintf("**********************************************\n"
            "Function:npt_event_solver_cop\n"
            "Estimate linear velocity for event cameras\n"
            "Author: Ji Zhao\n"
            "lastest update: 03.22.2025\n"
            "\n"
            "Input arguments\n"
            "1st: angular velocity (3*1 vector) or orientations (A two-level cell array, cell{cell{3x3 double}})\n"
            "2nd: events\n"
            "\n"
            "Output arguments\n"
            "1st: estimated linear velocity with scale ambiguity\n"
            "2nd: runtime (unit: microsecond)\n"
            "\n"
            "**********************************************\n"
            "\n");
        return;
    }

    // 1st input parameter
    // known angular velocity or orientations
    Vector3d omg;
    int numberOfVariables = 0;
    const mxArray *matlab_orientations;
    bool is_ang_vel_input = true;
    
    const mxArray *input1 = prhs[0];
    if (nrhs > 0 && !mxIsEmpty(input1)) {
        if (mxIsNumeric(input1) && (mxGetClassID(input1) == mxDOUBLE_CLASS)) {
            numberOfVariables = mxGetNumberOfElements(input1);
            if (numberOfVariables != 3) {
                mexErrMsgIdAndTxt("MATLAB:npt_event_solver_cop:invalidAngularVelocity", "First input vector must be 3x1.");
            }
            const double* omg0 = mxGetPr(input1);
            omg = Map<const Vector3d>(omg0);
        }
        else if (mxGetClassID(input1) == mxCELL_CLASS) {
            matlab_orientations = input1;
            is_ang_vel_input = false;
        }
    }
    else {
        mexErrMsgIdAndTxt("MATLAB:npt_event_solver_cop:invalidAngularVelocity", "First imput cannot be empty.");
    }

    // 2nd input parameter
    // events
    const mxArray *matlab_events;
    if (nrhs > 1 && !mxIsEmpty(prhs[1])) {
        matlab_events = prhs[1];
    }
    else {
        mexErrMsgIdAndTxt("MATLAB:npt_event_solver_cop:invalidEvents", "The second input cannot be empty.");
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

    // format conversion for orientations
    vector<vector<Matrix3d>> cpp_orientations;
    if (!is_ang_vel_input) {
        mwSize outerSize = mxGetNumberOfElements(matlab_orientations);
        for (mwSize i = 0; i < outerSize; ++i) {
            const mxArray* innerCell = mxGetCell(matlab_orientations, i);
            if (!mxIsCell(innerCell)) {
                mexErrMsgIdAndTxt("MATLAB:npt_event_solver_cop:invalidOrientations", "The first input of orintations is invalid.");
            }
            mwSize innerSize = mxGetNumberOfElements(innerCell);
            vector<Matrix3d> innerVec(innerSize);
            for (mwSize j = 0; j < innerSize; ++j) {
                const mxArray* matrixCell = mxGetCell(innerCell, j);
                innerVec[j] = Map<const Matrix<double, 3, 3, ColMajor>>(mxGetPr(matrixCell));
            }
            cpp_orientations.push_back(innerVec);
        }
    }

    // estimate translation
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
    Vector3d v = Vector3d::Zero();
    bool flag;
    if (is_ang_vel_input) {
        flag = compute_translation_cop(cpp_events, omg, v);
    }
    else {
        flag = translation_cop_core(cpp_events, cpp_orientations, v);
    }
    long long elapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
        double* result = mxGetPr(plhs[0]);
        for (int i = 0; i < 3; ++i) {
            result[i] = v(i);
        }
    }
    // Output the objective
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleScalar(static_cast<double>(elapsedTime));
    }

    return;
}
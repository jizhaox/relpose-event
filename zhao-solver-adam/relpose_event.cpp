#include <iostream>
#include <chrono>
#include "mex.h"
#include "matrix.h"
#include "optimizer_cop_flow_fdm.h"
#include "optimizer_cop_flow.h"
#include "optimizer_inc.h"
#include "optimizer_translation.h"

using namespace Eigen;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 0) {
        mexPrintf("**********************************************\n"
            "Function:relpose_event\n"
            "Estimate angular velocity and linear velocity for event cameras\n"
            "lastest update: 11.08.2024\n"
            "\n"
            "Input arguments\n"
            "1st: events\n"
            "2nd (optional): intial value for angular velocity\n"
            "3rd (optional): method types. Supported methods are listed below\n"
            "4th (optional): hyperparameters for formulations and optimiziers\n"
            "5th (optional): options for Adam optimzers\n"
            "\n"
            "Output arguments\n"
            "1st (optional): estimated angular velocity\n"
            "2nd (optional): estimated linear velocity with scale ambiguity\n"
            "3rd (optional): line structure (each cell correspnds to a line structure)\n"
            "4th (optional): objective\n"
            "5th (optional): runtime (unit: microsecond)\n"
            "\n"
            "Supported method types include:\n"
            "101: incidence + exact\n"
            "102: incidence + approximation\n"
            "103: incidence + cascade\n"
//            "201: CopRaw + exact + PGD\n"
//            "202: CopRaw + approximation + PGD\n"
//            "203: CopRaw + cascade + PGD\n"
//            "204: CopRaw + exact + AzEl\n"
//            "205: CopRaw + approximation + AzEl\n"
//            "206: CopRaw + cascade + AzEl\n"
            "301: coplanarity + exact\n"
            "302: coplanarity + approximation\n"
            "303: coplanarity + cascade\n"
            "304: coplanarity + exact + FDM\n"
            "305: coplanarity + approximation + FDM\n"
            "306: coplanarity + cascade + FDM\n"
            "\n"
            "FDM: finite difference method\n"
 //           "PGD: projected gradient descent\n"
 //           "AzEl: represent direction vector by azimuth and elevation angles\n"
            "**********************************************\n"
            "\n");
        return;
    }

    // 1st input parameter
    const mxArray *matlab_events;
    if (nrhs > 0 && !mxIsEmpty(prhs[0])) {
        matlab_events = prhs[0];
    }
    else {
        mexErrMsgIdAndTxt("MATLAB:relpose_event:invalidFirstInput", "The first input cannot be empty.");
    }
    
    // 2nd input parameter
    // intial value for the angular velocity
    Vector3d x = Vector3d::Zero();
    int numberOfVariables = 3;
    if (nrhs > 1 && !mxIsEmpty(prhs[1])) {
        numberOfVariables = mxGetNumberOfElements(prhs[1]);
        if (numberOfVariables != 3) {
            mexErrMsgIdAndTxt("MATLAB:relpose_event:invalidInitialValue", "Input vector must be 3x1.");
        }
        const double* omg0 = mxGetPr(prhs[1]);
        x = Map<const Vector3d>(omg0);
    }
    
    // 3rd input parameter
    // method type
    int method_type = 103;
    if (nrhs > 2 && !mxIsEmpty(prhs[2])) {
        method_type = static_cast<int>(mxGetScalar(prhs[2]));
    }

    // 4th input parameter
    // hyperparameters for problems and optimiziers
    double scale_lambda, diff_var, p_norm;
    scale_lambda = 1e6;
    diff_var = 1e-6; // Small deviation for numerical gradient computation
    p_norm = 1.0;
    if (nrhs > 3 && !mxIsEmpty(prhs[3])) {
        const mxArray* paramsStruct = prhs[3];
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

    // 5th input parameter
    // options for Adam optimzers
    double TolFun, TolX, MaxIter, MaxFunEvals, stepSize, beta1, beta2, epsilon, nEpochSize;
    TolFun = 1e-6;
    TolX = 1e-6;
    MaxIter = 1e6;
    MaxFunEvals = 1e6;
    stepSize = 0.001;
    beta1 = 0.9;
    beta2 = 0.999;
    epsilon = 1e-9;
    nEpochSize = 3;
    if (nrhs > 4 && !mxIsEmpty(prhs[4])) {
        const mxArray* optionsStruct = prhs[4];
        if (mxGetField(optionsStruct, 0, "TolFun") != NULL){
            TolFun = mxGetScalar(mxGetField(optionsStruct, 0, "TolFun"));
        }
        if (mxGetField(optionsStruct, 0, "TolX") != NULL){
            TolX = mxGetScalar(mxGetField(optionsStruct, 0, "TolX"));
        }
        if (mxGetField(optionsStruct, 0, "MaxIter") != NULL){
            MaxIter = static_cast<int>(mxGetScalar(mxGetField(optionsStruct, 0, "MaxIter")));
        }
        if (mxGetField(optionsStruct, 0, "MaxFunEvals") != NULL){
            MaxFunEvals = static_cast<int>(mxGetScalar(mxGetField(optionsStruct, 0, "MaxFunEvals")));
        }
        if (mxGetField(optionsStruct, 0, "stepSize") != NULL){
            stepSize = mxGetScalar(mxGetField(optionsStruct, 0, "stepSize"));
        }
        if (mxGetField(optionsStruct, 0, "beta1") != NULL){
            beta1 = mxGetScalar(mxGetField(optionsStruct, 0, "beta1"));
        }
        if (mxGetField(optionsStruct, 0, "beta2") != NULL){
            beta2 = mxGetScalar(mxGetField(optionsStruct, 0, "beta2"));
        }
        if (mxGetField(optionsStruct, 0, "epsilon") != NULL){
            epsilon = mxGetScalar(mxGetField(optionsStruct, 0, "epsilon"));
        }
        if (mxGetField(optionsStruct, 0, "nEpochSize") != NULL){
            nEpochSize = static_cast<int>(mxGetScalar(mxGetField(optionsStruct, 0, "nEpochSize")));
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

    // Do the job
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
    double fval = 0.0;
    if (method_type >= 101 && method_type <= 103)
    {
        fval = adam_inc(cpp_events, x, method_type, scale_lambda, diff_var, p_norm, 
            TolFun, TolX, MaxIter, MaxFunEvals, stepSize, beta1, beta2, epsilon, nEpochSize);
    }
//    else if (method_type >= 201 && method_type <= 203)
//    {
//        fval = adam_cop_raw_fdm(cpp_events, x, method_type, scale_lambda, diff_var, p_norm, 
//            TolFun, TolX, MaxIter, MaxFunEvals, stepSize, beta1, beta2, epsilon, nEpochSize);
//    }
//    else if (method_type >= 204 && method_type <= 206)
//    {
//        fval = adam_cop_raw_azel_fdm(cpp_events, x, method_type, scale_lambda, diff_var, p_norm, 
//            TolFun, TolX, MaxIter, MaxFunEvals, stepSize, beta1, beta2, epsilon, nEpochSize);
//    }
    else if (method_type >= 301 && method_type <= 303)
    {
        fval = adam_cop_flow(cpp_events, x, method_type, scale_lambda, diff_var, p_norm, 
            TolFun, TolX, MaxIter, MaxFunEvals, stepSize, beta1, beta2, epsilon, nEpochSize);
    }
    else if (method_type >= 304 && method_type <= 306)
    {
        fval = adam_cop_flow_fdm(cpp_events, x, method_type, scale_lambda, diff_var, p_norm, 
            TolFun, TolX, MaxIter, MaxFunEvals, stepSize, beta1, beta2, epsilon, nEpochSize);
    }
    else
    {
        mexErrMsgIdAndTxt("MATLAB:relpose_event:invalidMethodType", "Run relpose_event() to check supported methods.");
    }
    long long elapsedTime1 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count();

    // Output the optimized result
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(numberOfVariables, 1, mxREAL);
        double* result = mxGetPr(plhs[0]);
        for (int i = 0; i < numberOfVariables; ++i) {
            result[i] = x(i);
        }
    }

    // estimate translation and output
    long long elapsedTime2 = 0;
    vector<Matrix3d> line_struct_all;
    if (nlhs > 1) {
        startTime = std::chrono::high_resolution_clock::now();
        Vector3d v = Vector3d::Zero();
        bool flag = compute_translation_cop(cpp_events, x, v, line_struct_all);
        plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL);
        double* result2 = mxGetPr(plhs[1]);
        for (int i = 0; i < 3; ++i) {
            result2[i] = v(i);
        }
        elapsedTime2 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count();
    }

    // output line structure
    if (nlhs > 2) {
        size_t n = line_struct_all.size();
        plhs[2] = mxCreateCellMatrix(n, 1);
        for (size_t j = 0; j < n; ++j) {
            mxArray *mat = mxCreateDoubleMatrix(3, 3, mxREAL);
            double *data = mxGetPr(mat);
            const Eigen::Matrix3d& matrix = line_struct_all[j];
            std::memcpy(data, matrix.data(), 9 * sizeof(double));
            mxSetCell(plhs[2], j, mat);
        }
    }

    // output objective
    if (nlhs > 3) {
        plhs[3] = mxCreateDoubleScalar(fval);
    }

    // output runtime
    if (nlhs > 4) {
        // unit: microsecond
        plhs[4] = mxCreateDoubleMatrix(2, 1, mxREAL);
        double output3[2] = {static_cast<double>(elapsedTime1), static_cast<double>(elapsedTime2)};
        memcpy(mxGetData(plhs[4]), output3, 2 * sizeof(double));
    }
    return;
}
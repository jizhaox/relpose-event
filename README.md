## 1. Overview

This repository provides solvers to estimate egomotion for event cameras. It contains solvers for angular and linear velocities  [1]. The core optimization algorithms are implemented in C++. We provide compiled MEX files compatible with both Ubuntu 20.04 + Matlab 2018b, as well as Windows 10 + Matlab R2022a. The data generation and demonstration scripts are written in MATLAB.


```
@inproceedings{zhao2025full,
	title={Full-DoF Egomotion Estimation for Event Cameras Using Geometric Solvers},
	author={Zhao, Ji and Guan, Banglei and Liu, Zibin and Kneip, Laurent},
	booktitle={IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
	year={2025}
}
```

**Author:** [Ji Zhao](https://sites.google.com/site/drjizhao) and [Zibin Liu](https://github.com/Zibin6)

## 2. Quick Start
 
run `test_2dof_solver.m` to estimate linear velocity only.

run `test_5dof_solver_adam.m` to esimtate angular and linear velocities.

If you want to compile the MEX file, change the path for EIGEN library in file `compile.m` and run it.
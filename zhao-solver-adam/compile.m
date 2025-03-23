if isunix
    % Linux or macOS
    eigenPath = '/usr/include/eigen3';
    mex(['-I' eigenPath], 'CXXFLAGS=-std=c++14 -fPIC', '-O', 'relpose_event.cpp');
elseif ispc
    % Windows
    eigenPath = 'D:/eigen-3.4.0';
    mex(['-I' eigenPath], 'CXXFLAGS=$CXXFLAGS /std:c++14', ...
        'CXXOPTIMFLAGS=$CXXOPTIMFLAGS /O2 /DNDEBUG', 'relpose_event.cpp');
    mex(['-I' eigenPath], 'CXXFLAGS=$CXXFLAGS /std:c++14', ...
        'CXXOPTIMFLAGS=$CXXOPTIMFLAGS /O2 /DNDEBUG', 'npt_event_solver_cop.cpp');
else
    disp('Unknown operating system');
end

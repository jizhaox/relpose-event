if isunix
    % Linux
    mex('-I/usr/include/eigen3', 'CXXFLAGS=-std=c++14 -fPIC', '-O', 'relpose_event.cpp');
elseif ispc
    % Windows
    mex -ID:/BaiduNetdiskDownload/eigen-3.4.0 CXXFLAGS='$CXXFLAGS /std:c++14' CXXOPTIMFLAGS='\$CXXOPTIMFLAGS /O2 /DNDEBUG' relpose_event.cpp
else
    disp('Unknown operating system');
end

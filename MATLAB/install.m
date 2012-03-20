% This is a small script to compile the necessary MEX files.
cflags  = '-O3 -Os -std=gnu99';
ldflags = '-O3 -Os -s';
options = sprintf('-O -largeArrayDims COPTIMFLAGS="%s" LDOPTIMFLAGS="%s"',...
		  cflags,ldflags);
sourcefiles = {' ../C/doublevectormatlab.c'
	       ' ../C/singlematrixmatlab.c' 
	       ' ../C/sigmoid.c'
	       ' ../C/vectorops.c'
	       ' ../C/diagsq.c'
	       ' ../C/diagsqt.c'
	       ' ../C/varbvs.c'
	       ' ../C/varbvsbin.c' };
eval(['mex ',options,' ../C/var1matlab.c',sourcefiles{[1 2 4]}]);
eval(['mex ',options,' ../C/varbvsupdatematlab.c',sourcefiles{[1:4 7]}]);
eval(['mex ',options,' ../C/diagsqmatlab.c',sourcefiles{[1 2 4 5]}]);
fprintf('Compilation of MEX files is complete.\n');

% mex -O -largeArrayDims CXXOPTIMFLAGS='-O3' LDOPTIMFLAGS='-O3' ...
%     diagsqtfast.cpp doublevector.cpp singlematrix.cpp common.cpp
% mex -O -largeArrayDims CXXOPTIMFLAGS='-O3' LDOPTIMFLAGS='-O3' ...
%     varbvsbinupdate.cpp doublevector.cpp singlematrix.cpp common.cpp

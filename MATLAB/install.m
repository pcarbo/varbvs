% This is a small script to compile the necessary MEX files.

% Here is the opportunity to set some additional flags here that may be
% passed to the compiler. These flags tell the GCC compiler to use the ISO
% C99 standard, and to optimize the code as much as possible. Depending
% on the compiler you use to build the MEX shared library files, you may
% want to change these flags, or set them to the empty string ('').
cflags  = '-std=gnu99 -O3 -Os';
ldflags = '-s -O3 -Os';

% These are the files containing the main functions implemented in C. Note
% that not all these files are needed to compile each of the MEX files.
corefiles = {' ../C/doublevectormatlab.c'
	       ' ../C/singlematrixmatlab.c' 
	       ' ../C/sigmoid.c'
	       ' ../C/vectorops.c'
	       ' ../C/diagsq.c'
	       ' ../C/diagsqt.c'
	       ' ../C/varbvs.c'
	       ' ../C/varbvsbin.c' };

% These are the commands to build the build the MEX shared library files.
options = sprintf('-O -largeArrayDims COPTIMFLAGS="%s" LDOPTIMFLAGS="%s"',...
		  cflags,ldflags);
eval(['mex ',options,' ../C/var1matlab.c',corefiles{[1 2 4]}]);
return
eval(['mex ',options,' ../C/varbvsupdatematlab.c',corefiles{[1:4 7]}]);
eval(['mex ',options,' ../C/diagsqmatlab.c',corefiles{[1 2 4 5]}]);
fprintf('Compilation of MEX files is complete.\n');

% mex -O -largeArrayDims CXXOPTIMFLAGS='-O3' LDOPTIMFLAGS='-O3' ...
%     diagsqtfast.cpp doublevector.cpp singlematrix.cpp common.cpp
% mex -O -largeArrayDims CXXOPTIMFLAGS='-O3' LDOPTIMFLAGS='-O3' ...
%     varbvsbinupdate.cpp doublevector.cpp singlematrix.cpp common.cpp

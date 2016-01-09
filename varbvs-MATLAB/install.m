% This script compiles the necessary MEX files.

% These are the files containing the main functions implemented in C.
corefiles = { 'doublevectormex.c '
              'doublematrixmex.c '
	      'singlematrixmex.c ' 
	      '../varbvs-R/src/vectorops.c '
	      '../varbvs-R/src/sigmoid.c '
	      '../varbvs-R/src/varbvs.c '
	      '../varbvs-R/src/varbvsbin.c ' };

% These are the commands to build the build the MEX shared library files.
opts = sprintf('-O -largeArrayDims -DMATLAB_MEX_FILE -I../varbvs-R/src ');
eval(['mex ',opts,'diagsqmex.c ',corefiles{1:4}]);
eval(['mex ',opts,'diagsqtmex.c ',corefiles{1:4}]);
eval(['mex ',opts,'var1mex.c ',corefiles{1:4}]);
eval(['mex ',opts,'varbvsnormupdatemex.c ',corefiles{1:6}]);
eval(['mex ',opts,'varbvsbinupdatemex.c ',corefiles{[1:5 7]}]);
eval(['mex ',opts,'varbvsbinzupdatemex.c ',corefiles{[1:5 7]}]);
fprintf('Compilation of MEX files is complete.\n');

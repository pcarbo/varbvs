% This script compiles the necessary MEX files.

% These are the commands to build the build the MEX shared library files.
opts = sprintf('-O -largeArrayDims -DMATLAB_MEX_FILE -I../varbvs-R/src');

% Build diagsqmex MEX file.
% eval(['mex ',opts,'diagsqmex.c doublevectormex.c singlematrixmex.c ',...
%       '../varbvs-R/src/misc.c');

% Build diagsqtmex MEX file.
eval(['mex ',opts,' diagsqtmex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c']);

% Build var1mex MEX file.
eval(['mex ',opts,' var1mex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c']);

% Build varbvsnormupdatemex MEX file.
eval(['mex ',opts,' varbvsnormupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c ../varbvs-R/src/misc.c ../varbvs-R/src/varbvs.c']);

% Build varbvsbinupdatemex MEX file.
eval(['mex ',opts,' varbvsbinupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c ../varbvs-R/src/misc.c ../varbvs-R/src/varbvsbin.c']);

% Build varbvsbinzupdatemex MEX file.
eval(['mex ',opts,' varbvsbinzupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c doublematrixmex.c ../varbvs-R/src/misc.c ',...
      '../varbvs-R/src/varbvsbin.c']);
fprintf('Compilation of MEX files is complete.\n');

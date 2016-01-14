% This script compiles the necessary MEX files.

% These are the commands to build the build the MEX shared library files.
opts = sprintf('-O -largeArrayDims -DMATLAB_MEX_FILE -I../varbvs-R/src');

% Build diagsqmex MEX file.
fprintf('Building diagsqmex MEX file.\n');
eval(['mex ',opts,' diagsqmex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c ../varbvs-R/src/diagsq.c']);
fprintf('\n');

% Build diagsqtmex MEX file.
fprintf('Building diagsqtmex MEX file.\n');
eval(['mex ',opts,' diagsqtmex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c ../varbvs-R/src/diagsq.c']);
fprintf('\n');

% Build var1mex MEX file.
fprintf('Building var1mex MEX file.\n');
eval(['mex ',opts,' var1mex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c']);
fprintf('\n');

% Build varbvsnormupdatemex MEX file.
fprintf('Building varbvsnormupdatemex MEX file.\n');
eval(['mex ',opts,' varbvsnormupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c ../varbvs-R/src/misc.c ../varbvs-R/src/varbvs.c']);
fprintf('\n');

% Build varbvsbinupdatemex MEX file.
fprintf('Building varbvsbinupdatemex MEX file.\n');
eval(['mex ',opts,' varbvsbinupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c ../varbvs-R/src/misc.c ../varbvs-R/src/varbvs.c']);
fprintf('\n');

% Build varbvsbinzupdatemex MEX file.
fprintf('Building varbvsbinzupdatemex MEX file.\n');
eval(['mex ',opts,' varbvsbinzupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c doublematrixmex.c ../varbvs-R/src/misc.c ',...
      '../varbvs-R/src/varbvs.c']);
fprintf('\n');

fprintf('Compilation of MEX files is complete.\n');

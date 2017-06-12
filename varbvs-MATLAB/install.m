% This script compiles the necessary MEX files.

% These are the commands to build the build the MEX shared library files.
% I've found that the -std=gnu99 flag is useful for successfully building
% the MEX files in MATLAB 8.1 (R2013a), but may not be needed for other
% MATLAB releases.

% Part of the varbvs package, https://github.com/pcarbo/varbvs
%
% Copyright (C) 2012-2017, Peter Carbonetto
%
% This program is free software: you can redistribute it under the
% terms of the GNU General Public License; either version 3 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANY; without even the implied warranty of
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
%
opts = ['-O -largeArrayDims -DMATLAB_MEX_FILE -I../varbvs-R/src ' ...
        'COPTIMFLAGS="-std=gnu99"'];

% Build diagsqmex MEX file.
fprintf('Building diagsqmex MEX file.\n');
eval(['mex ',opts,' diagsqmex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c ../varbvs-R/src/diagsq.c']);

% Build diagsqtmex MEX file.
fprintf('Building diagsqtmex MEX file.\n');
eval(['mex ',opts,' diagsqtmex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c ../varbvs-R/src/diagsq.c']);

% Build var1mex MEX file.
fprintf('Building var1mex MEX file.\n');
eval(['mex ',opts,' var1mex.c doublevectormex.c singlematrixmex.c ',...
      '../varbvs-R/src/misc.c']);

% Build varbvsnormupdatemex MEX file.
fprintf('Building varbvsnormupdatemex MEX file.\n');
eval(['mex ',opts,' varbvsnormupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c ../varbvs-R/src/misc.c ../varbvs-R/src/varbvs.c']);

% Build varbvsbinupdatemex MEX file.
fprintf('Building varbvsbinupdatemex MEX file.\n');
eval(['mex ',opts,' varbvsbinupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c ../varbvs-R/src/misc.c ../varbvs-R/src/varbvs.c']);

% Build varbvsbinzupdatemex MEX file.
fprintf('Building varbvsbinzupdatemex MEX file.\n');
eval(['mex ',opts,' varbvsbinzupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c doublematrixmex.c ../varbvs-R/src/misc.c ',...
      '../varbvs-R/src/varbvs.c']);

% Build varbvsmixupdatemex MEX file.
fprintf('Building varbvsmixupdatemex MEX file.\n');
eval(['mex ',opts,' varbvsmixupdatemex.c doublevectormex.c ',...
      'singlematrixmex.c doublematrixmex.c ../varbvs-R/src/misc.c ',...
      '../varbvs-R/src/varbvs.c']);

fprintf('Compilation of MEX files is complete.\n');

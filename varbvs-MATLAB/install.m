% This is a small script to compile the necessary MEX files.

% Here is the opportunity to set some additional flags here that may be
% passed to the compiler. These flags tell the GCC compiler to use the ISO
% C99 standard, and to optimize the code as much as possible. Depending
% on the compiler you use to build the MEX shared library files, you may
% want to change these variables, or set them to the empty string ('').
cflags  = '-std=gnu99 -O3 -Os';
ldflags = '-O3 -Os';

% These are the files containing the main functions implemented in C. Note
% that not all these files are needed to compile each of the MEX files.
Rsrcdir   = '../R/varbvs/src/';
corefiles = {'C/doublevectormex.c '
             'C/doublematrixmex.c '
	     'C/singlematrixmex.c ' 
	     [ Rsrcdir 'vectorops.c ' ]
	     [ Rsrcdir 'sigmoid.c '   ]
	     [ Rsrcdir 'varbvs.c '    ]
             [ Rsrcdir 'varbvsmix.c ' ]
	     [ Rsrcdir 'varbvsbin.c ' ]};

% These are the commands to build the build the MEX shared library files.
options = sprintf(['-O -largeArrayDims -IC -I%s ' ...
		   'COPTIMFLAGS="%s" LDOPTIMFLAGS="%s" '],...
		   Rsrcdir,cflags,ldflags);
eval(['mex ',options,'C/var1mex.c ',corefiles{1:4}]);
eval(['mex ',options,'C/diagsqmex.c ',corefiles{1:4}]);
eval(['mex ',options,'C/diagsqtmex.c ',corefiles{1:4}]);
eval(['mex ',options,'C/varbvsupdatemex.c ',corefiles{1:6}]);
eval(['mex ',options,'C/varbvsmixupdatemex.c ',corefiles{[1:5 7]}]);
eval(['mex ',options,'C/varbvsbinupdatemex.c ',corefiles{[1:5 8]}]);
eval(['mex ',options,'C/varbvszbinupdatemex.c ',corefiles{[1:5 8]}]);
fprintf('Compilation of MEX files is complete.\n');

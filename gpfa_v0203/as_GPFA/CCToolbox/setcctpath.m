function SetCCPath()

base = fileparts(mfilename('fullpath'));
addpath(genpath(base), '-begin');
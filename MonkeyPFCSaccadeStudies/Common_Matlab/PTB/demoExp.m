% I am trying to figure out how to get Adam's experiments up and running.

proj_path = fullfile('/Users', 'chadwickboulay', 'Desktop', 'SachsLab', 'MatlabProjects');
addpath(proj_path);
addpath(fullfile(proj_path, 'fp_files'));
addpath(fullfile(proj_path, 'PsychoPhysics_DAQ_commands'));
addpath(fullfile(proj_path, 'Radial Spatial task'));
addpath(fullfile(proj_path, 'Spatial task'));
addpath(fullfile(proj_path, 'online calibration'));
addpath(genpath(fullfile(proj_path, 'Saccade GUI')));
clear proj_path

Saccade_GUI;

% 

% calls saccadeGUImain;


%(p)ause, z-new trial block, (q)uit
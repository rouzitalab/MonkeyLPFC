function showmodel(model,trajs);
%SHOWMODEL  Plot learned cluster model parameters and/or data
%  SHOWMODEL(model,[trajs]) will plot the learned model parameters in 'model'
%  using VIEWMODEL. If 'trajs' exists, then this curve data will be plotted also.
%
%  This function should be considered a hybrid function/script. This file is
%  meant to be directly edited by the user to enable various kinds of
%  plots. Feel free to change all of the 'ops' lines to suit your needs.

% Scott J Gaffney   28 February 2004
% Computer Science
% University of California, Irvine

PROGNAME = 'showmodel';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Style options that don't need to be changed
%
% In general you don't want to change these
ops.PlotLines              = 1;   % turn on-off line plotting
ops.PlotColors             = 1;   % turn on-off color plotting
ops.PlotSymbols            = 0;   % turn on-off symbol plotting

% Specifiy a default style used for all classes and cyclones.
% Empty signifies to plot class-specific styles.
ops.linespec               = '-';
ops.colorspec              = [];
ops.symbolspec             = [];

% Change the order of styles that are used
ops.Colors                 = 'rgbkmc';
ops.Styles                 = {'-','-.','--',':'};
ops.Symbols                = 'xosdv^<>*+ph';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General plotting options that are usually changed.
%
ops.PlotMean               = 1;   % plot mean curves
ops.plot_regcoef.len       = 15;  % length of mean curves
ops.LineWidth              = 2;   % width of mean curves
ops.MeanCurveColor         = '';  % plot means as white, for example
ops.MarkerSize             = 8;   % change marker size
ops.FontSize               = 12;  % change axes font size

% Special options for plotting cyclones (4 measurements per day)
ops.PlotMeansDayMarkers    = 0;   % plot day markers on mean curves
ops.PlotTrajsDayMarkers    = 0;   % plot day markers on data curves
ops.PlotMarkerFaceColor    = 0;   % fill in day markers

% Selective curve plotting
ops.SelectClass            = [];  % select specific classes to plot
ops.SeparateFigs           = 0;   % plot clusters in separate figures
ops.SelectTrajs            = [];  % select specific curves to plot
ops.PlotNumRand            = [];  % plot random number of curves



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIEWMODEL specific options
%
ops.ViewModel              = 1;   % enable plotting of models
ops.PlotAlign              = 1;   % plot using alignments

% Selecting types of plots (you can plot all of them at once)
ops.TwoD                   = 0;   % e.g., plot y_1 vs. y_2
ops.PlotDimensions         = 1;   % plot y_1 vs. time,  y_2 vs. time, etc.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the plotting

% view the curve models in relative space
if (ops.ViewModel)
  viewmodel(model, trajs, ops);
  vfig = gcf;
end









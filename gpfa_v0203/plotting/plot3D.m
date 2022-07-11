function plot3D(seq, target, xspec, varargin)
%
% plot3D(seq, xspec, ...)
%
% Plot neural trajectories in a three-dimensional space.
%
% INPUTS:
%
% seq        - data structure containing extracted trajectories
% xspec      - field name of trajectories in 'seq' to be plotted 
%              (e.g., 'xorth' or 'xsm')
%
% OPTIONAL ARGUMENTS:
%
% dimsToPlot - selects three dimensions in seq.(xspec) to plot 
%              (default: 1:3)
% nPlotMax   - maximum number of trials to plot (default: 20)
% redTrials  - vector of trialIds whose trajectories are plotted in red
%              (default: [])
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  dimsToPlot = 1:3;
  nPlotMax   = 20;
  redTrials  = [];
  assignopts(who, varargin);

  if size(seq(1).(xspec), 1) < 3
    fprintf('ERROR: Trajectories have less than 3 dimensions.\n');
    return
  end

  f = figure;
  pos = get(gcf, 'position');
%   set(f, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
  col = [[1 0 0]; [0 1 0]; [0 0 1]; [0 1 1]; [1 1 0]; [1 0 1];
      [0.5 0.1 0.1]; [0.1 0.5 0.1]; [0.5 0.5 0.5]; [0.5 0.5 1];
      [1 0.5 0.5]; [0.3 0.2 0.4]];
  for n = 1:min(length(seq), nPlotMax)
    dat = seq(n).(xspec)(dimsToPlot,:);
    T   = seq(n).T;
        
%     if ismember(seq(n).trialId, redTrials)
%       col = [1 0 0]; % red
%       lw  = 3;
%     else
%       col = 0.8 * [1 1 1]; % gray
%       lw = 0.1;
%     end
    plot3(dat(1,:), dat(2,:), dat(3,:), '-', 'color', col(target(n)+1,:));
    hold on;
%     plot3(dat(1,6), dat(2,6), dat(3,6), 'o', 'linewidth', 2.5, 'color', [1 0 0]);
%     hold on;
%     plot3(dat(1,9), dat(2,9), dat(3,9), 'o', 'linewidth', 2.5, 'color', [0 1 0]);
%     hold on;
%     plot3(dat(1,19), dat(2,19), dat(3,19), 'o', 'linewidth', 2.5, 'color', [0 0 1]);
%     hold on;
  end

  axis equal;
  if isequal(xspec, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(2));
    str3 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  else
    str1 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(2));
    str3 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  end
  xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
  ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);
  zlabel(str3, 'interpreter', 'latex', 'fontsize', 24);

function plot2D(seq, y, xspec,titl, varargin)
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
%              (default: 1:2)
% nPlotMax   - maximum number of trials to plot (default: 20)
% redTrials  - vector of trialIds whose trajectories are plotted in red
%              (default: [])
%
% @ 2009 Byron Yu -- byronyu@stanford.edu
  dimsToPlot = [1 2];
  nPlotMax   = 400;
  redTrials  = [];
  assignopts(who, varargin);

%   if size(seq(1).(xspec), 1) < 2
%     fprintf('ERROR: Trajectories have less than 2 dimensions.\n');
%     return
%   end
  
  figure;
%   pos = get(gcf, 'position');
%   set(f, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
%   col = [[1 0 0]; [0 1 0]; [0 0 1]; [0 1 1]; [1 1 0]; [1 0 1];
%       [0.5 0.1 0.1]; [0.1 0.5 0.1]; [0.5 0.5 0.5]; [0.5 0.5 1];
%       [1 0.5 0.5]; [0.3 0.2 0.4]];
  for n = 1:min(length(seq), nPlotMax)
    dat = seq(n).(xspec)(dimsToPlot,:);
%     dat = tmp - tmp(:,1);
    if y(n) == 'r'
        col = [1 0 0];
    elseif y(n) == 'g'
        col = [0 1 0];
    else
        col = [0 0 1];
    end
    T   = seq(n).T;
        
%     if ismember(seq(n).trialId, redTrials)
%       col = [1 0 0]; % red
%       lw  = 3;
%     else
%       col = 0.8 * [1 1 1]; % gray
%       lw = 0.1;
%     end
%     plot(dat(1,:),'.-','color', col(y(n)+1,:));
%     hold on;
    plot(dat(1,:), dat(2,:), '.-','color', col);
    hold on;
%     plot(dat(1,30), dat(2,30), 'o', 'color', col(target(n)+1,:));
%     hold on;
%     plot(dat(1,43), dat(2,43), 'o', 'color', col(target(n)+1,:));
%     hold on;
%     plot(dat(1,93), dat(2,93), 'o', 'color', col(target(n)+1,:));
%     hold on;
%     plot(dat(1,end), dat(2,end), 'o', 'color', col(y(n)+1,:));
%     hold on;
%     plot(dat(1,end), dat(2,end), 'o', 'linewidth', 1.5, 'color', col);
%     hold on;
%     plot(dat(1,9), dat(2,9), 'o', 'linewidth', .1, 'color', [0.5 0.5 0.5]);
%     hold on;
%     plot(dat(1,19), dat(2,19), 'o', 'linewidth', .1, 'color', [0.5 0.5 0.5]);
  end

  axis equal;
  if isequal(xspec, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(2));
%     str3 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  else
    str1 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(2));
%     str3 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  end
  xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
  ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);
  title(titl);
%   zlabel(str3, 'interpreter', 'latex', 'fontsize', 24);

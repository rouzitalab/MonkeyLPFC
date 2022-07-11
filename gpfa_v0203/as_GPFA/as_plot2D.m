function outIdx=as_plot2D(trialDirLookup, seq, xspec, varargin)
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

  dimsToPlot = 1:2;
  nPlotMax   = min(80,size(seq,2)); % ajs July 2011 20;
  redTrials  = [];
  assignopts(who, varargin);

  if size(seq(1).(xspec), 1) < 3
    fprintf('ERROR: Trajectories have less than 3 dimensions.\n');
    return
  end

  f = figure;
%   pos = get(gcf, 'position');
%   set(f, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
  
  outIdx=[];
  for n = 1:min(length(seq), nPlotMax)
    dat = seq(n).(xspec)(dimsToPlot,:);
    T   = seq(n).T;
        
    if ismember(seq(n).trialId, redTrials)
      col = [1 0 0]; % red
      lw  = 3;
    elseif trialDirLookup(seq(n).trialId,2)== 1
      col =  [1 0 0]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    elseif trialDirLookup(seq(n).trialId,2)== 2
      col =  [0 1 0]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    elseif trialDirLookup(seq(n).trialId,2)== 3
      col =  [0 0 1]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    elseif trialDirLookup(seq(n).trialId,2)== 4
      col = [0 1 1]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    elseif trialDirLookup(seq(n).trialId,2)== 5
      col =  [1 0 1]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    elseif trialDirLookup(seq(n).trialId,2)== 6
      col =  [1 1 0]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    elseif trialDirLookup(seq(n).trialId,2)== 7
      col =  [1 1 1]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    elseif trialDirLookup(seq(n).trialId,2)== 8
      col =  [.2 .4 1]; % gray
      lw = 0.5;
      outIdx(n)=trialDirLookup(seq(n).trialId,1);
    end
    plot(dat(1,:), dat(2,:), '.-', 'linewidth', lw, 'color', col);
    hold on;
  end

  axis equal;
  if isequal(xspec, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(2));
  else
    str1 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(2));
  end
  xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
  ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);


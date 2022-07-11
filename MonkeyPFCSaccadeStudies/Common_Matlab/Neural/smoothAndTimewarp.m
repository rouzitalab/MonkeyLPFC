function [trajectories, traj_timestamp, win_durs] = smoothAndTimewarp(cDat, analysisParams, tw_wins, dhParams, dhRedMode, winSamples)
%[trajectories, traj_timestamp, win_durs] = smoothAndTimewarp(cDat, analysisParams, tw_wins, dhParams, dhRedMode, winSamples)

nNeurons = sum(dhParams.keep_neurons);

%% Extract smoothed (maybe dim reduced) trajectories from each window
cDat = trimTrials(cDat, analysisParams.anaWins);
[~, long_tr_ix] = max(arrayfun(@(x)(sum(x.timeBool)), cDat.trial));
full_tvec = cDat.trial(long_tr_ix).tVec(cDat.trial(long_tr_ix).timeBool);

nRed = nNeurons;
if dhRedMode ~= -1
    nRed = 8;
end
D = cDat2DataHigh(cDat.trial);
[reducedD, ~, ~, ~] =...
        reducedims(D, dhRedMode, nRed, dhParams);

%% Calculate the common longest timestamp vector for the reduced data.
% We are going to use the middle of each bin as its timestamp.
reduced_nsamps = arrayfun(@(x)(size(x.data, 2)), reducedD);
% reduced_tvec = dhParams.step_size*(0:max(reduced_nsamps)-1);
reduced_tvec = full_tvec(1:dhParams.step_size:end);
reduced_tvec = reduced_tvec(1:max(reduced_nsamps))...
    + ceil(dhParams.binWidth / 2);

%% Collect window fragments, time-warping when necessary.
nWins = length(tw_wins);
red_out = cell(length(reducedD), nWins);
out_tvec = cell(1, nWins);
win_durs = nan(1, nWins);

% test_ind = zeros(length(cDat.trial), length(reduced_tvec));

for win_ix = 1:nWins
    ana_win = tw_wins(win_ix);
    
    % Calculate the in-window time stamps for each trial (saved into
    % .trial.tVec(.timeBool))
    win_cdat = trimTrials(cDat, ana_win);
    
    % Get the reducedD data that overlaps with this window
    win_nsamples = nan(size(reducedD));
    for tr_ix = 1:length(reducedD)
        this_trial = win_cdat.trial(tr_ix);
        this_tvec = cDat.trial(tr_ix).tVec(this_trial.timeBool);
        this_rtvec = reduced_tvec(1:reduced_nsamps(tr_ix));
        this_rbool = this_rtvec >= this_tvec(1) & this_rtvec <= this_tvec(end);
        red_out{tr_ix, win_ix} = reducedD(tr_ix).data(:, this_rbool);
        win_nsamples(tr_ix) = sum(this_rbool);
        
%         test_ind(tr_ix, find(this_rbool)) = test_ind(tr_ix, find(this_rbool)) + win_ix;
    end
    
    % Calculate the target number of samples to time-warp to
    if exist('winSamples', 'var') && length(winSamples)==nWins
        targetSamples = winSamples(win_ix);
    else
        targetSamples = ceil(median(win_nsamples));
    end
    win_durs(win_ix) = targetSamples;
    out_tvec{win_ix} = ana_win.plotX(1) + dhParams.step_size*(0:targetSamples-1);
    
    % Time-warp to targetSamples
    uq_nsamples = unique(win_nsamples);
    if length(uq_nsamples) > 1  % any(win_nsamples ~= targetSamples)
        % Prepare time-vectors for each possible number of samples.
        uq_tvec = cell(1, length(uq_nsamples));
        for uq_ix = 1:length(uq_nsamples)
            uq_tvec{uq_ix} = 1:(targetSamples-1)/(uq_nsamples(uq_ix)-1):targetSamples;
        end
        
        % Linear interpolation of windows to stretch or compress so
        % they're all the same size.
        for tr_ix = 1:length(reducedD)
            if win_nsamples(tr_ix) ~= targetSamples
                tr_tvec = uq_tvec{uq_nsamples==win_nsamples(tr_ix)};
                % Unfortunately linterp only operates on vectors.
                % Do each channel separately.
                outD = nan(nRed, targetSamples);
                for chan_ix = 1:nRed
                    outD(chan_ix, :) = linterp(tr_tvec, red_out{tr_ix, win_ix}(chan_ix, :), 1:targetSamples);
                end
                red_out{tr_ix, win_ix} = outD;
                clear tr_tvec outD chan_ix
            end
        end
    end
    
end
% imagesc(test_ind)
%% Return the trajectories in one big matrix.
trajectories = nan(size(red_out, 1), sum(win_durs), nRed);
for tr_ix = 1:size(red_out, 1)
    trajectories(tr_ix, :, :) = cat(2, red_out{tr_ix, :})';
end
traj_timestamp = cat(2, out_tvec{:});
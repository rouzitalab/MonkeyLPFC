% Single Unit Activities
% base = 0;
% for ch = 1:32
% if ch < 10
% wf = load(dPath + "waveforms_00" + ch + ".mat").waveforms;
% else
% wf = load(dPath + "waveforms_0" + ch + ".mat").waveforms;
% end
% w_idx = sort(nonzeros(sorts(ch).maxRatings));
% w_idx = [0; w_idx];
% for trial = 1:nTrials
% for unit = 1:length(w_idx)
% idx = find((wf.trialInfo.trial==trial)&(wf.units==w_idx(unit)));
% for spike = 1:length(idx)
% sua{trial, base+unit} = [sua{trial,base+unit}, wf.waves(w_idx(unit)+1, idx(spike))];
% end
% end
% end
% base = base + length(w_idx);
% end
%%%%%%%%%%% WORK WITH TIMES %%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
nUnits = 0;
dPath = "G:\Projects\MonkeyLPFC\sorted_analyses\mk\JerryLee\090623\sra3_1_j_050_00+\";
% dPath = "/users/alireza/projects/MonkeyLPFC/Data/sra3_1_j_050_00+/";
sorts = load(dPath + "sorts.mat").sorts;
for ch = 1:32
    if ch < 10
        wf = load(dPath + "waveforms_00" + ch + ".mat").waveforms;
    else
        wf = load(dPath + "waveforms_0" + ch + ".mat").waveforms;
    end
    w_idx = sort(nonzeros(sorts(ch).maxRatings));
    w_idx = [0; w_idx];
    nu = length(w_idx);
    nUnits = nUnits+nu;
end
conditions = wf.trialInfo.condition;
nTrials = length(conditions);
s = wf.trialInfo.trialStartTimes;
e = wf.trialInfo.trialEndTimes;
i = [];
j = [];
v = [];
base = 1;
for ch = 1:32
    if ch < 10
        wf = load(dPath + "waveforms_00" + ch + ".mat").waveforms;
    else
        wf = load(dPath + "waveforms_0" + ch + ".mat").waveforms;
    end
    w_idx = sort(nonzeros(sorts(ch).maxRatings));
    w_idx = [0; w_idx];
    nu = length(w_idx);
    for samp = 1:length(wf.spikeTimes)
        if ~isnan(wf.trialInfo.trial(samp))
            i(end+1) = base+wf.units(samp);
            j(end+1) = ceil(wf.spikeTimes(samp)-s(1));
            v(end+1) = 1;
        else
            continue
        end
    end
    base = base + nu;
end
sua = sparse(i,j,v);

clc
clear
dPath = "G:\Projects\MonkeyLPFC\sorted_analyses\mk\JerryLee\090623\sra3_1_j_050_00+\";
% dPath = "/users/alireza/projects/MonkeyLPFC/Data/sra3_1_j_050_00+/";
sorts = load(dPath + "sorts.mat").sorts;
units = zeros(1,32);
for ch = 1:32
    if ch < 10
        wf = load(dPath + "waveforms_00" + ch + ".mat").waveforms;
    else
        wf = load(dPath + "waveforms_0" + ch + ".mat").waveforms;
    end
    w_idx = sort(nonzeros(sorts(ch).maxRatings));
    w_idx = [0; w_idx];
    units(ch) = length(w_idx);
end
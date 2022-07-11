clear;
clc;
path = 'G:\Projects\neuroGAM-master\data';
data_list = dir(fullfile(path, '*.mat'));
filename = 'session_analysis.xlsx';
col_name = 'CDEFGHIJKLMNOPQRSTUVWXYZ';
for sess = 2 : 9
clearvars -except 'sess' 'data_list' 'path' 'filename' 'col_name'
fprintf('Working on Sess %d\n', sess);
data_path = append(path, '\',data_list(sess).name);
data1 = load(data_path);
learned = data1.learned + 1;
unlearned = data1.unlearned + 1;
target = data1.target;
cue = data1.color;
outcome = data1.outcome;
color = double.empty();
color(length(cue)) = 0;
for i=1:length(cue) % assigning labels to the cue color : Red = 1, Green = 2, Blue = 3
    if isequal(cue{i},'r')
        color(i)=1;
    elseif isequal(cue{i},'g')
        color(i)=2;
    elseif isequal(cue{i},'b')
        color(i)=3;
    end
end
trials = length(target); % number of time trials
rul = zeros(1, trials);
for i = 1:trials
    if target(i)==0 && color(i)==1, rul(i)=1;
    elseif target(i)==0 && color(i)==2, rul(i)=2;
    elseif target(i)==0 && color(i)==3, rul(i)=3;
    elseif target(i)==1 && color(i)==1, rul(i)=4;
    elseif target(i)==1 && color(i)==2, rul(i)=5;
    elseif target(i)==1 && color(i)==3, rul(i)=6;
    elseif target(i)==2 && color(i)==1, rul(i)=7;
    elseif target(i)==2 && color(i)==2, rul(i)=8;
    elseif target(i)==2 && color(i)==3, rul(i)=9;
    elseif target(i)==3 && color(i)==1, rul(i)=10;
    elseif target(i)==3 && color(i)==2, rul(i)=11;
    elseif target(i)==3 && color(i)==3, rul(i)=12;
    elseif target(i)==4 && color(i)==1, rul(i)=13;
    elseif target(i)==4 && color(i)==2, rul(i)=14;
    elseif target(i)==4 && color(i)==3, rul(i)=15;
    elseif target(i)==5 && color(i)==1, rul(i)=16;
    elseif target(i)==5 && color(i)==2, rul(i)=17;
    elseif target(i)==5 && color(i)==3, rul(i)=18;
    elseif target(i)==6 && color(i)==1, rul(i)=19;
    elseif target(i)==6 && color(i)==2, rul(i)=20;
    elseif target(i)==6 && color(i)==3, rul(i)=21;
    elseif target(i)==7 && color(i)==1, rul(i)=22;
    elseif target(i)==7 && color(i)==2, rul(i)=23;
    else, rul(i)=24;
    end
end
cor = find(outcome==0);
icor = find(outcome==9);
c_l=[];
ic_ul=[];
for i = 1 : length(learned)
    if sum(cor==learned(i))
        c_l = [c_l, learned(i)];
    end
end
for i = 1 : length(unlearned)
    if sum(icor==unlearned(i))
        ic_ul = [ic_ul, unlearned(i)];
    end
end
for i = 1 : 24
    cel = [col_name(i), int2str(4*sess-6)];
    writematrix(sum(rul(learned)==i),filename,'Sheet',1,'Range',cel);
    cel = [col_name(i), int2str(4*sess-5)];
    writematrix(sum(rul(c_l)==i),filename,'Sheet',1,'Range',cel);
    cel = [col_name(i), int2str(4*sess-4)];
    writematrix(sum(rul(unlearned)==i),filename,'Sheet',1,'Range',cel);
    cel = [col_name(i), int2str(4*sess-3)];
    writematrix(sum(rul(ic_ul)==i),filename,'Sheet',1,'Range',cel);
end
fprintf('Session %d Done!\n\n', sess);
end
fprintf('All Done!\n');
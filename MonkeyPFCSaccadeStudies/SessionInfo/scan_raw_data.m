%This script examines the data files to see which files match up with each
%other.
addpath(genpath(fullfile(pwd, '..', 'Common')));
my_paths; global paths
my_subjects = my_subs;

%% Get info of nev files
output = struct('subname', '', 'dirname', '', 'fname', '',...
    'fileID', '', 'trialInfo', '', 'ns2', false, 'ccf', false, 'datenum', []);
row_ix = 0;
for sub_ix = 1:length(my_subjects)
    this_sub = my_subjects(sub_ix);
    dirsub = dir(fullfile(paths.cdataRoot, this_sub.name));
    dirsub = dirsub([dirsub.isdir]);
    dirsub = dirsub(~ismember({dirsub.name}, {'.', '..'}));
    for dir_ix = 1:length(dirsub)
        dirnev = dir(fullfile(paths.cdataRoot, this_sub.name, dirsub(dir_ix).name, '*.nev'));
        for file_ix = 1:length(dirnev)
            row_ix = row_ix + 1;
            fname = fullfile(paths.cdataRoot, this_sub.name, dirsub(dir_ix).name, dirnev(file_ix).name);
            fileInfo = getNevFileInfo(fname);
            output(row_ix).subname = this_sub.name;
            output(row_ix).dirname = dirsub(dir_ix).name;
            output(row_ix).fname = dirnev(file_ix).name;
            output(row_ix).fileID = fileInfo.fileID;
            output(row_ix).trialInfo = fileInfo.trialInfo;
            output(row_ix).nTrials = length(fileInfo.trialInfo);
            output(row_ix).ns2 = exist(fullfile(paths.cdataRoot, this_sub.name, dirsub(dir_ix).name, [dirnev(file_ix).name(1:end-3) 'ns2']), 'file') == 2;
            output(row_ix).ccf = exist(fullfile(paths.cdataRoot, this_sub.name, dirsub(dir_ix).name, [dirnev(file_ix).name(1:end-3) 'ccf']), 'file') == 2;
            output(row_ix).datenum = dirnev(file_ix).datenum;
        end
    end
end
nevFileInfo = output;
clear row_ix sub_ix this_sub dirsub dir_ix dirnev file_ix fname fileInfo output

%% Get info of ptb files
output = struct('nCentreOut', [], 'nCueColours', [], 'nDiagonals', [],...
    'subname', '', 'fname', '', 'fileID', [], 'datenum', [],...
    'nTrials', [], 'eDataInfo', [],...
    'pFilename', '', 'pTrials', [], 'pExpCode', [], 'pExpInfo', [],...
    'edfFile', '');
row_ix = 0;
for sub_ix = 1:length(my_subjects)
    this_sub = my_subjects(sub_ix);
    dirptb = dir(fullfile(paths.ptbedfRoot, this_sub.sname, 'ptbmat', '*.ptbmat'));
    for file_ix = 1:length(dirptb)
        row_ix = row_ix + 1;
        ptb = getPTB(fullfile(paths.ptbedfRoot, this_sub.sname, 'ptbmat', dirptb(file_ix).name));
        %TODO: Get experiment type info from ptb file.
        exp_type = nan(length(ptb.eData.trial), 1);
        for t_ix = 1:length(ptb.eData.trial)
            tr = ptb.eData.trial(t_ix);
            
            if tr.expType == 'C' || tr.expType == 'D'
                % Training-level. See
                % FromAdam/SachsLab/MatlabProjects/Saccade
                % GUI/SR3/dsr3T2c_setscreens.m around line 187.
                if isfield(tr, 'SR3trainingLevel') && ~isempty(tr.SR3trainingLevel)
                    trl = tr.SR3trainingLevel;
                elseif isfield(tr, 'SR3InitialtrainingLevel') && ~isempty(tr.SR3InitialtrainingLevel)
                    trl = tr.SR3InitialtrainingLevel;
                else
                    trl = ptb.params.SR3trainingLevel;
                end
                
                if isfield(tr, 'SR3errorStrategy')  && ~isempty(tr.SR3errorStrategy)
                    es = tr.SR3errorStrategy;
                elseif isfield(ptb.params, 'SR3errorStrategy')
                    es = ptb.params.SR3errorStrategy;
                else
                    es = nan;
                end
                
                if tr.expType == 'D' && length(trl) > 1
                    trl = trl(2);
                    es = es(2);
                end
                
                if isstruct(es) 
                    if strcmpi(es.error, 'resample')
                        es = 7;
                    elseif strcmpi(es.error, 'repeat')
                        es = 8;
                    elseif strcmpi(es.error, 'probation')
                        es = 9;
                    end
                end
                
                if trl == 3 || trl == 4
                    exp_type(t_ix) = 1;  % Centre-out
                elseif tr.expType == 'C' && trl == 0
                    exp_type(t_ix) = 2;  % Cue-colour centre-out
                elseif tr.expType == 'D' && trl == 0 && es == 9
                    exp_type(t_ix) = 2;  % Cue-colour centre-out
                elseif tr.expType == 'D' && trl == 0 && es ~= 9
                    exp_type(t_ix) = 3;  % Diagonals;
                end
            end
        end
        
        output(row_ix).nCentreOut = sum(exp_type == 1);
        output(row_ix).nCueColours = sum(exp_type == 2);
        output(row_ix).nDiagonals = sum(exp_type == 3);
            
        
        output(row_ix).subname = this_sub.name;
        output(row_ix).fname = dirptb(file_ix).name;
        output(row_ix).fileID = dirptb(file_ix).name(1:end-7);
        output(row_ix).datenum = dirptb(file_ix).datenum;
        
        eDataInfo = ptb.eData;
        output(row_ix).nTrials = length(eDataInfo.trial);
        eDataInfo = rmfield(eDataInfo, 'trial');
        output(row_ix).eDataInfo = eDataInfo;
        
        
        output(row_ix).pFilename = ptb.params.filename;
        output(row_ix).pTrials = ptb.params.numberofTrials;
        output(row_ix).pExpCode = ptb.params.expCode;
        output(row_ix).pExpInfo = ptb.params.expInfo;
        
        diredf = dir(fullfile(paths.ptbedfRoot, this_sub.sname, 'edf', [output(row_ix).fileID '*.edf']));
        if ~isempty(diredf)
            output(row_ix).edfFile = diredf.name;
        end
        
        fprintf('Found %i/%i trials in %s dated %s.\n',...
            output(row_ix).nTrials, output(row_ix).pTrials,...
            output(row_ix).fname,...
            datestr(output(row_ix).datenum));     
    end
end
ptbFileInfo = output;
clear output row_ix sub_ix this_sub dirptb file_ix ptb exp_type t_ix
clear tr trl es exp_type eDataInfo diredf

%% Try to merge info - unique fileIDs for both nev and ptb that match.

ptbIndForNev = nan(length(nevFileInfo), 1);

ptbIDs = {ptbFileInfo.fileID};
nevIDs = {nevFileInfo.fileID};

ptbIsUnique = false(length(ptbIDs), 1);
for pp = 1:length(ptbIDs)
    ptbIsUnique(pp) = sum(strcmpi(ptbIDs, ptbIDs{pp})) == 1;
end

for nn = 1:length(nevIDs)
    if sum(strcmpi(nevIDs, nevIDs{nn})) == 1
        if any(strcmpi(ptbIDs(ptbIsUnique), nevIDs{nn}))
            ptbIndForNev(nn) = find(strcmpi(ptbIDs, nevIDs{nn}));
        end
    end
end
clear pp nn ptbIsUnique

%Check
% nbool = ~isnan(ptbIndForNev);
% pbool = ptbIndForNev(~isnan(ptbIndForNev));
% [nevIDs(nbool)' {nevFileInfo(nbool).nTrials}'...
%     ptbIDs(pbool)' {ptbFileInfo(pbool).nTrials}'...
%     num2cell([nevFileInfo(nbool).nTrials] - [ptbFileInfo(pbool).nTrials])']
%I found one that was off. sra3_2_j_047_00+02... but the nev file seems
%large enough. Maybe it was just a digital communication error.
% nidx = find(nbool);
% badidx = nidx(abs([nevFileInfo(nbool).nTrials] - [ptbFileInfo(pbool).nTrials]) > 1);
% ptbIndForNev(badidx) = nan;
% clear nbool pbool badidx nidx

%% Merge info - those that FileID did not match, try date and trial number match
% leftOvers = nevFileInfo(isnan(ptbIndForNev));
% [{leftOvers.subname}' {leftOvers.dirname}' {leftOvers.fname}' {leftOvers.fileID}' {leftOvers.nTrials}']
maxtdiff = 10/(24*60); %10 minutes
for nn = 1:length(ptbIndForNev)
    if isnan(ptbIndForNev(nn))
        nfi = nevFileInfo(nn);
        pbool = abs([ptbFileInfo.nTrials] - nfi.nTrials) <= 1;
        pbool = pbool & abs([ptbFileInfo.datenum] - nfi.datenum) <= maxtdiff;
        if sum(pbool)==1
            ptbIndForNev(nn) = find(pbool);
        end
    end
end
clear maxtdiff nn nfi pbool
%This leaves me without matches for sum(isnan(ptbIndForNev)) == 116
%Of those 31 files with more than 100 trials...
%find(isnan(ptbIndForNev) & [nevFileInfo.nTrials]' >= 100)'
%Behavioural data are missing from June 2 - 10 (nn = 145:164) then neural
%data missing from June 10 to June 17 (but there are behavioural data from
%June 11,15,16)
%sra3_2_j_0filename_00+01.ptbmat = 2000 trials, April 16
%sra3_2_j_017_00+.mat = 1880 trials, April 21 ... but I'm missing nev files
%from that date. A folder is there...
%sra3_2_j_017_00+.ptbmat = 1561 trials, April 23; no nev files.

%% Build session table
header = {'Subject' 'Date' 'ExpCode' 'IsGood' 'NTrials' 'NEVDir', 'NEVFName' 'NS2FName' 'PTBFName' 'EDFName' 'NCentreOut' 'NCueColour' 'NDiagonal'};
nind = find(~isnan(ptbIndForNev))';
output = cell(length(nind), length(header));
for nn = 1:length(nind)
%     nevFileInfo(nind(nn))
%     ptbFileInfo(ptbIndForNev(nind(nn)))
    nfi = nevFileInfo(nind(nn));
    pfi = ptbFileInfo(ptbIndForNev(nind(nn)));
    output{nn, strcmpi(header, 'Subject')} = nfi.subname;
    output{nn, strcmpi(header, 'Date')} = datestr(nfi.datenum, 29);%
    output{nn, strcmpi(header, 'ExpCode')} = pfi.pExpCode;
    output{nn, strcmpi(header, 'IsGood')} = 1;
    output{nn, strcmpi(header, 'NTrials')} = nfi.nTrials;
    output{nn, strcmpi(header, 'NEVDir')} = nfi.dirname;
    output{nn, strcmpi(header, 'NEVFName')} = nfi.fname;
    if nfi.ns2
        output{nn, strcmpi(header, 'NS2FName')} = [nfi.fname(1:end-3) 'ns2'];
    else
        output{nn, strcmpi(header, 'NS2FName')} = 'none';
    end
    %output{nn, strcmpi(header, 'HasCCF')} = nfi.ccf;
    %output{nn, strcmpi(header, 'PTBFileID')} = char(pfi.eDataInfo.fileID);
    output{nn, strcmpi(header, 'PTBFName')} = pfi.fname;
    output{nn, strcmpi(header, 'EDFName')} = pfi.edfFile;
    output{nn, strcmpi(header, 'NCentreOut')} = pfi.nCentreOut;
    output{nn, strcmpi(header, 'NCueColour')} = pfi.nCueColours;
    output{nn, strcmpi(header, 'NDiagonal')} = pfi.nDiagonals;
end
clear nind nn nfi pfi
%%
if ~isdir(paths.sessionInfo)
    mkdir(paths.sessionInfo);
end
fname = fullfile(paths.sessionInfo, 'SessionList.csv');
fid = fopen(fname, 'w');
fprintf(fid, [strjoin(header,',') '\n']);
for row_ix = 1:size(output,1)
    fprintf(fid, '%s,%s,%s,%i,%i,%s,%s,%s,%s,%s,%i,%i,%i\n', output{row_ix,:});
end
fclose(fid);
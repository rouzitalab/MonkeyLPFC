%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
sessions = my_sessions('Contains_Inverse');
subjects = my_subs();
global paths

%% Quick init for debug.
sess_ix = find(strcmpi({sessions.ptbfname}, 'sra3_1_j_061_00+02.ptbmat'));
this_sess = sessions(sess_ix);

%% Continue pre-processing --> save cDat
for sess_ix = 1:length(sessions)
    this_sess = sessions(sess_ix);
    
    %% TODO: Figure out how big backPorch can be
    
    % Load data.
    cDat = getCerbDat(this_sess, 'neur_pref', {'cDat', 'mksort'});
    
    % save cDat
    if ~isdir(fullfile(paths.preprocessed, 'cDat'))
        mkdir(fullfile(paths.preprocessed, 'cDat'));
    end
    saveFull = fullfile(paths.preprocessed, 'cDat', [this_sess.edfname(1:end-3) 'mat']);
    save(saveFull, '-struct', 'cDat');
end
clear sess_ix this_sess cDat saveFull
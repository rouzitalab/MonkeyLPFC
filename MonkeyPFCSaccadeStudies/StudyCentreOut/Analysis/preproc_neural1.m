
%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
sessions = my_sessions('CentreOut');
subjects = my_subs();
global paths

addpath(genpath(fullfile(paths.ml3rd, 'mksort')));


%% Extract all the data to mksort format.
for sess_ix = 1:length(sessions)
    this_sess = sessions(sess_ix);

    %%
    nevdir = fullfile(paths.cdataRoot, this_sess.subject, this_sess.nevdir);
    mksortdir = fullfile(paths.preprocessed, 'mk', this_sess.subject, ...
        this_sess.nevdir, this_sess.ptbfname(1:end-7));
    if ~isdir(mksortdir)
        mkdir(mksortdir);
    end

    %%
    %Transfer the data from nev to mksort
    nev2MatWaveforms(nevdir, mksortdir,{{this_sess.nevfname}}); %, artifactThresh, snippetLen, firstSnippetPt, conserveMemory)

    %% Add trial/condition info
    %Trial information will not be used for sorting, but it will be useful to
    %get it into the mksort output, because then we can ignore the nev file and
    %skip the nevmat file.
    mergeNevPtbIntoMKSort(fullfile(nevdir, this_sess.nevfname),...
        fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname),...
        mksortdir);

end
clear sess_ix this_sess nevdir
%% DO Manual sorting
%Use mksort
%Load the created directory.
%Sort.
%Save.
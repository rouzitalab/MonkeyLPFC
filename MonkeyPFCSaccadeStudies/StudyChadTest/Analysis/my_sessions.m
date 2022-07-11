function sessions = my_sessions(varargin)
% sessions = my_sessions(fav_name)
% 
% fav_name is an optional string argument to limit the sessions returned to
% those in a set of favourites defined below.
%
% sessions is an array of structs containing fields that match the columns
% in SessionInfo/SessionList.csv

% Below you should define a block of sessions you wish to analyze.
% Each block is a structure with fields
% name      - A convenient name to describe the sessions to be included.
% sessions  - A cell array of string triplets. The first string is the name
%             that matches the Subject column, the second matches NEVDir,
%             and the third matches PTBFname.

favourites.name = 'TrickySessions';  % Multiple experiment types per file.
favourites.sessions = [...
    {'Marty' '090814' 'sra3_1_m_065_00+03'};...  %
    {'Marty' '090818' 'sra3_1_m_067_00+01'};...  %A
    {'Marty' '090819' 'sra3_1_m_068_00+01'};...  %B
    {'Marty' '090902' 'sra3_1_m_069_00+01'};...  %C
    {'Marty' '090908' 'sra3_1_m_070_00+02'}];    %A

sessions = getSessions();

%% If we were provided with a list of favourites, reduce to those
if nargin > 0 && sum(strcmpi({favourites.name}, varargin{1})) == 1
    fav_sess = favourites(strcmpi({favourites.name}, varargin{1})).sessions;
    sess_bool = false(size(sessions));
    for sess_ix = 1:size(fav_sess, 1)
        sess_bool = sess_bool | ...
            (strcmpi({sessions.subject}', fav_sess{sess_ix, 1}) &...
            strcmpi({sessions.nevdir}', fav_sess{sess_ix, 2}) &...
            strcmpi({sessions.ptbfname}', [fav_sess{sess_ix, 3} '.ptbmat']));
    end
    sessions = sessions(sess_bool);
end


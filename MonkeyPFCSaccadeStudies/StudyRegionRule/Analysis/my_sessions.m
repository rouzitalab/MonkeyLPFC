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

favourites.name = 'RegionRule';
favourites.sessions = [...
    {'JerryLee' '090630' 'sra3_2_j_054_00+02'};...  %B; many blocks, good performance. s33
    {'JerryLee' '090702' 'sra3_2_j_056_00+02'};...  %A; many blocks, good performance. s21
    {'JerryLee' '090703' 'sra3_2_j_057_00+01'};...  %B; many blocks, good performance. s35
    {'JerryLee' '090708' 'sra3_2_j_060_00+01'};...  %B; many blocks, good performance. s23
    {'JerryLee' '090716' 'sra3_2_j_066_00+01'}];%...  %B; many blocks, good performance. s30
    
    
% Good but not used
%     {'JerryLee' '090626' 'sra3_2_j_053_00+01'};...  %A; blue ur, ok performance. s27
%     {'JerryLee' '090701' 'sra3_2_j_055_00+01'};...  %C; many blocks, good performance
%     {'JerryLee' '090706' 'sra3_2_j_057_00+06'};...  %C; many blocks, good performance
%     {'JerryLee' '090707' 'sra3_2_j_059_00+02'};...  %A; many blocks, good performance
%     {'JerryLee' '090709' 'sra3_2_j_061_00+01'};...  %C; many blocks, good performance
%     {'JerryLee' '090710' 'sra3_2_j_062_00+02'};...  %A; many blocks, good performance
%     {'JerryLee' '090711' 'sra3_2_j_063_00+01'};...  %B; OK, not great performance
%     {'JerryLee' '090712' 'sra3_2_j_064_00+01'};...  %C; many blocks, good performance
%     {'JerryLee' '090715' 'sra3_2_j_065_00+01'};...  %A; many blocks, good performance
%     {'JerryLee' '090717' 'sra3_2_j_067_00+01'};...  %C; many blocks, good performance
%     {'Marty' '091012' 'sra3_2_m_093_00+02'};...     %A; red ur, blue dl
%     {'Marty' '091027' 'sra3_2_m_100_00+01'};...     %C; r dr
%     {'Marty' '091030' 'sra3_2_m_101_00+05'};...     %B; b dl
%     {'Marty' '091106' 'sra3_2_m_105_00+01'}];       %C; g dr, b ul
    


% Bad
%     {'Marty' '091004' 'sra3_2_m_087_00+03'};...     %B; poor performance, combines with next
%     {'Marty' '091004' 'sra3_2_m_087_00+04'};...     %B; poor performance, for red ul
%     {'Marty' '091005' 'sra3_2_m_088_00+02'};...     %C; no overlap
%     {'Marty' '091009' 'sra3_2_m_090_00+01'};...     %A; no overlap
%     {'Marty' '091011' 'sra3_2_m_092_00+04'};...     %B; only perfect overlap
%     {'Marty' '091016' 'sra3_2_m_094_00+04'};...     %C; green dl, blue ur, but performance bad
%     {'Marty' '091017' 'sra3_2_m_095_00+01'};...     %C; no overlap
%     {'Marty' '091018' 'sra3_2_m_096_00+01'};...     %C; no overlap
%     {'Marty' '091025' 'sra3_2_m_098_00+01'};...     %A; not much overlap r dl, g ul
%     {'Marty' '091026' 'sra3_2_m_099_00+01'};...     %B; only perfect overlap blue up up up
%     {'Marty' '091103' 'sra3_2_m_104_00+01'};...     %B; no overlap
%     {'Marty' '091107' 'sra3_2_m_106_00+01'}];       %B; not much overlap

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


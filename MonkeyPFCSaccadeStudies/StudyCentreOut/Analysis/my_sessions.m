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

favourites.name = 'CentreOut';
favourites.sessions =[...    
    {'JerryLee' '090617' 'sra3_2_j_047_00+01'};...	%A:21,431; 80/85*
    {'JerryLee' '090626' 'sra3_1_j_053_00+'};...    %A:26,229; 89/79*
    {'JerryLee' '090701' 'sra3_1_j_054_00+'};...    %C:17,259; 64/48*
    {'JerryLee' '090709' 'sra3_1_j_061_00+02'};...  %C:14,234; 61/53*
    {'JerryLee' '090711' 'sra3_1_j_063_00+02'};...  %B:27,258; 84/89*
    {'JerryLee' '090716' 'sra3_1_j_066_00+01'};...  %B:23,266; 87/71*
    {'Marty' '090708' 'sra3_1_m_035_00+03'};...     %B:14,231; 56/30*
    {'Marty' '090712' 'sra3_1_m_039_00+01'};...     %C:15,184; 52/47*
    {'Marty' '090713' 'sra3_1_m_040_00+02'};...     %A:21,228; 58/46*
    {'Marty' '090716' 'sra3_1_m_043_00+01'};...     %A:28,328; 62/55*
    {'Marty' '090723' 'sra3_1_m_050_00+02'};...     %B:13,200; 44/51
    {'Marty' '090724' 'sra3_1_m_051_00+01'}];       %C:20,200; 59/50*
    
% Good but not used
% {'JerryLee' '090618' 'sra3_2_j_048_00+01'};...  %B:26,300; 62/92
% {'JerryLee' '090702' 'sra3_1_j_056_00+'};...    %A:18,273; 73/77
% {'JerryLee' '090708' 'sra3_1_j_060_00+01'};...  %B:15,245; 77/80
% {'JerryLee' '090712' 'sra3_1_j_064_00+01'};...  %A:18,206; 57/55
% {'JerryLee' '090715' 'sra3_1_j_065_00+01'};...  %A:23,254; 76/89
% {'Marty' '090710' 'sra3_1_m_037_00+01'};...       %A:22,257; 55/44
% {'Marty' '090711' 'sra3_1_m_038_00+02'};...     %B:16-2,239; 46/77
% {'Marty' '090717' 'sra3_1_m_044_00+01'};...     %B:17,186; 41/38*
% {'Marty' '090718' 'sra3_1_m_045_00+03'};...     %C:21,190; 35/40
% {'Marty' '090723' 'sra3_1_m_050_00+02'};...     %B:13,200; 44/51
% {'Marty' '090726' 'sra3_1_m_053_00+01'};...     %B:17,191; 35/47
% {'Marty' '090730' 'sra3_1_m_056_00+01'};...     %B:8,189; 42/39
% {'Marty' '091006' 'sra3_1_m_089_00+01'};...       %A:14,405; 57/45
% {'Marty' '091009' 'sra3_1_m_090_00+02'};...       %A:13,580; 50/39
% {'Marty' '091012' 'sra3_1_m_093_00+01'};...       %A:18,484; 52/32. eye tracker drift. Maybe ok.
% {'Marty' '091018' 'sra3_1_m_096_00+01'}];%;...  %C:10,453; 23/29
% {'Marty' '090719' 'sra3_1_m_046_00+01'};...       %A:?,296; doesn't exist in csv.

% Bad
% {'JerryLee' '090622' 'sra3_2_j_049_00+01'};...  %C:24,245. Bad behavioural data for target 1.
% {'JerryLee' '090717' 'sra3_1_j_067_00+01'}];%;...  %C:?,192; bad behav
% {'JerryLee' '090707' 'sra3_1_j_059_00+01'};...  %A:?,350; bad behav
% {'JerryLee' '090710' 'sra3_1_j_062_00+01'};...  %A:?,302; bad behav
% {'Marty' '090714' 'sra3_1_m_041_00+01'};...       %B:~7,304; difficult to sort, a few really noisy events
% {'Marty' '090725' 'sra3_1_m_052_00+03'};...       %A:?,203; bad behaviour
% {'Marty' '090731' 'sra3_1_m_057_00+02'};...       %A:?,198; bad behaviour (only UR/DL)
% {'Marty' '090911' 'sra3_1_m_072_00+02'};...       %B:?,260; bad behaviour (only UR/DL; not many)
% {'Marty' '091004' 'sra3_1_m_087_00+01'};...       %B:?noisy,367; eye tracker drift. Maybe ok.
% {'Marty' '091016' 'sra3_1_m_094_00+01'};...  %C:9,502; bad behav 9/502. Only l-r and ur-dl
% {'Marty' '091017' 'sra3_1_m_095_00+01'};...  %C:?,531 pretty noisy
% {'Marty' '091025' 'sra3_1_m_098_00+01'};...       %A:?noisy units,452
% {'Marty' '091026' 'sra3_1_m_099_00+01'};...       %B:,418; bad behaviour. Very few UL/DR
% {'Marty' '091027' 'sra3_1_m_100_00+01'};...  %C:?,402 pretty noisy
% {'Marty' '091011' 'sra3_1_m_092_00+01'};...       %B:,638; behaviour empty?
% {'Marty' '091030' 'sra3_1_m_101_00+01'};...       %B:,504; bad behaviour? (empty)
% {'Marty' '091103' 'sra3_1_m_104_00+01'};...       %B:,401; units not great
% {'Marty' '091106' 'sra3_1_m_105_00+01'};...       %C:,394; very flat units
% {'Marty' '091107' 'sra3_1_m_106_00+01'}];%;...       %B:,501; units not great

% sessions = getSessions();
sessions = ["sra3_2_j_047_00+01","sra3_1_j_053_00+","sra3_1_j_054_00+","sra3_1_j_061_00+02","sra3_1_j_063_00+02","sra3_1_j_066_00+01","sra3_1_m_035_00+03","sra3_1_m_039_00+01","sra3_1_m_040_00+02","sra3_1_m_043_00+01","sra3_1_m_050_00+02","sra3_1_m_051_00+01"];
%% If we were provided with a list of favourites, reduce to those
% if nargin > 0 && sum(strcmpi({favourites.name}, varargin{1})) == 1
%     fav_sess = favourites(strcmpi({favourites.name}, varargin{1})).sessions;
%     sess_bool = false(size(sessions));
%     for sess_ix = 1:size(fav_sess, 1)
%         sess_bool = sess_bool | ...
%             (strcmpi({sessions.subject}', fav_sess{sess_ix, 1}) &...
%             strcmpi({sessions.nevdir}', fav_sess{sess_ix, 2}) &...
%             strcmpi({sessions.ptbfname}', [fav_sess{sess_ix, 3} '.ptbmat']));
%     end
%     sessions = sessions(sess_bool);
% end


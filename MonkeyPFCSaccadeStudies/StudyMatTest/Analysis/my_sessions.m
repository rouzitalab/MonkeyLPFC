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
favourites.name = 'Contains_Inverse';
favourites.sessions = [...
    {'JerryLee' '090622' 'sra3_2_j_049_00+02'};...  %C; s26. OK! 1 complement.
    {'JerryLee' '090623' 'sra3_1_j_050_00+'};...    %A; s33. Good performance. 8 blocks. Complements.
    {'JerryLee' '090624' 'sra3_1_j_051_00+'}];%...    %B; Good! 8 blocks. Complements.
% 
% favourites(end+1).name = 'Blocks';
% favourites(end).sessions = [...
%     {'JerryLee' '090521' 'sra3_2_j_035_00+02'};...  %C; Performance waivers. No complements
%     {'JerryLee' '090522' 'sra3_2_j_036_00+01'};...  %C-part1. OK performance. No complements.
%     {'JerryLee' '090522' 'sra3_2_j_036_00+02'};...  %C-part2. OK. 1 complement.
%     {'JerryLee' '090601' 'sra3_2_j_037_00+03'};...  %A; OK. 1 complement.
%     {'JerryLee' '090617' 'sra3_2_j_047_00+02'};...  %A; OK. 1 complement.
%     {'JerryLee' '090618' 'sra3_2_j_048_00+02'};...  %B; no complements. Performance at end is bad.
%     {'JerryLee' '090622' 'sra3_2_j_049_00+02'};...  %C; OK! 1 complement.
%     {'JerryLee' '090623' 'sra3_1_j_050_00+'};...    %A; Good! 8 blocks. Complements.
%     {'JerryLee' '090624' 'sra3_1_j_051_00+'};...    %B; Good! 8 blocks. Complements.
%     {'Marty' '090707' 'sra3_1_m_034_00+03'};...  %A - only one block
%     {'Marty' '090708' 'sra3_1_m_035_00+04'};...  %B - only one block
%     {'Marty' '090710' 'sra3_1_m_037_00+02'};...  %A - only one block
%     {'Marty' '090711' 'sra3_1_m_038_00+03'};...  %B - only one block
%     {'Marty' '090712' 'sra3_1_m_039_00+02'};...  %C - only one block
%     {'Marty' '090713' 'sra3_1_m_040_00+03'};...  %A - only one block
%     {'Marty' '090714' 'sra3_1_m_041_00+03'};...  %B - only one block
%     {'Marty' '090715' 'sra3_1_m_042_00+04'};...  %C - only one block
%     {'Marty' '090716' 'sra3_1_m_043_00+02'};...  %A - only one block
%     {'Marty' '090717' 'sra3_1_m_044_00+02'};...  %B - only one block
%     {'Marty' '090718' 'sra3_1_m_045_00+04'};...  %C - only one block
%     {'Marty' '090719' 'sra3_1_m_046_00+02'};...  %A. NevDir out of order. - only one block
%     {'Marty' '090721' 'sra3_1_m_048_00+04'};...  %C - only one block
%     {'Marty' '090722' 'sra3_1_m_049_00+02'};...  %A - only one block
%     {'Marty' '090723' 'sra3_1_m_050_00+03'};...  %B - only one block
%     {'Marty' '090724' 'sra3_1_m_051_00+02'};...  %C - only one block
%     {'Marty' '090725' 'sra3_1_m_052_00+04'};...  %A - only one block
%     {'Marty' '090726' 'sra3_1_m_053_00+02'};...  %B - only one block
%     {'Marty' '090730' 'sra3_1_m_056_00+02'};...  %B - only one block
%     {'Marty' '090731' 'sra3_1_m_057_00+02'};...  %A - only one block
%     {'Marty' '090806' 'sra3_1_m_058_00+01'}];%   %C - only one block
% 
% %Listed sessions are complete for JL
% %Marty has many more sessions, but files are mixed experiment types.
% favourites(end+1).name = 'All';
% favourites(end).sessions = [...
%     {'JerryLee' '090427' 'sra3_2_j_018_00+'};...    %A - Only 2 classes and I don't get its rule.
%     {'JerryLee' '090521' 'sra3_2_j_035_00+02'};...  %C; Performance waivers. No complements
%     {'JerryLee' '090522' 'sra3_2_j_036_00+01'};...  %C-part1. OK performance. No complements.
%     {'JerryLee' '090522' 'sra3_2_j_036_00+02'};...  %C-part2. s5. OK. 1 complement. bad unit data.
%     {'JerryLee' '090601' 'sra3_2_j_037_00+03'};...  %A; s10. OK. 1 complement. bad unit data.
%     {'JerryLee' '090617' 'sra3_2_j_047_00+02'};...  %A; s24. OK. 1 complement. Only 167 trials.
%     {'JerryLee' '090618' 'sra3_2_j_048_00+02'};...  %B; no complements. Performance at end is bad.
%     {'JerryLee' '090622' 'sra3_2_j_049_00+02'};...  %C; OK! 1 complement.
%     {'JerryLee' '090623' 'sra3_1_j_050_00+'};...    %A; Good! 8 blocks. Complements.
%     {'JerryLee' '090624' 'sra3_1_j_051_00+'};...    %B; Good! 8 blocks. Complements.
%     {'Marty' '090707' 'sra3_1_m_034_00+03'};...  %A - only one block
%     {'Marty' '090708' 'sra3_1_m_035_00+04'};...  %B - only one block
%     {'Marty' '090710' 'sra3_1_m_037_00+02'};...  %A - only one block
%     {'Marty' '090711' 'sra3_1_m_038_00+03'};...  %B - only one block
%     {'Marty' '090712' 'sra3_1_m_039_00+02'};...  %C - only one block
%     {'Marty' '090713' 'sra3_1_m_040_00+03'};...  %A - only one block
%     {'Marty' '090714' 'sra3_1_m_041_00+03'};...  %B - only one block
%     {'Marty' '090715' 'sra3_1_m_042_00+04'};...  %C - only one block
%     {'Marty' '090716' 'sra3_1_m_043_00+02'};...  %A - only one block
%     {'Marty' '090717' 'sra3_1_m_044_00+02'};...  %B - only one block
%     {'Marty' '090718' 'sra3_1_m_045_00+04'};...  %C - only one block
%     {'Marty' '090719' 'sra3_1_m_046_00+02'};...  %A. NevDir out of order. - only one block
%     {'Marty' '090721' 'sra3_1_m_048_00+04'};...  %C - only one block
%     {'Marty' '090722' 'sra3_1_m_049_00+02'};...  %A - only one block
%     {'Marty' '090723' 'sra3_1_m_050_00+03'};...  %B - only one block
%     {'Marty' '090724' 'sra3_1_m_051_00+02'};...  %C - only one block
%     {'Marty' '090725' 'sra3_1_m_052_00+04'};...  %A - only one block
%     {'Marty' '090726' 'sra3_1_m_053_00+02'};...  %B - only one block
%     {'Marty' '090730' 'sra3_1_m_056_00+02'};...  %B - only one block
%     {'Marty' '090731' 'sra3_1_m_057_00+02'};...  %A - only one block
%     {'Marty' '090806' 'sra3_1_m_058_00+01'}];%   %C - only one block

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


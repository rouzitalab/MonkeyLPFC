function [S FILENAME PATHNAME] = loadnev(FILENAME)
% Load spikes from nev file
%   S = loadnev will open a dialog to chose the nev file and return a
%   structure with the spikes.
%
%   S = loadnev(FILENAME) will skip the open dialog window.
%
%   The result S is an array of structure, each element containing:
%       S{i}.Channel: The channel identifier of this element
%       
%       S{i}.Spikes : A NxM matrices containing N spikes of M samples
%
%       S{i}.Times: A Nx1 vector containing the time at which the spikes
%       occured. The time is in samples (independent of sampling frequency)
%
%   The code is optimized for speed, not memory usage. Thus, it shouldn't
%   be used on a 32bits machine nor a machine with less than 6GB of memory.
%   Of course this depends on the actual size of the file to load.
%
%   Authors: Hicham Semmaoui <hicham.semmaoui@polymtl.ca>
%            Jonathan Drolet <jonathan.drolet@polymtl.ca>
%            Polystim Laboratory www.polystim.ca
%
%   Original authors: Unknown
if nargout < 1
    error('No output no work!');
end;

if nargin == 0
    [FILENAME,PATHNAME]=uigetfile({'*.nev'}, 'Choose nev file to convert');
    if FILENAME==0
        error('No filename given');
    end;
elseif nargin == 1
    PATHNAME = '';
end;

%Waitbar ------------------------------------------------------------------
Wb = waitbar (0, 'Please wait, loading file...');
%--------------------------------------------------------------------------
fid=fopen([PATHNAME FILENAME],'r');

if fid == 0
    error('Invalid filename.');
end;

%read general header
FileType=fread(fid,8,'*char');
Version=fread(fid,2,'uchar');
FileFormatAdditional=fread(fid,2,'*char');
HeaderSize=fread(fid,1,'uint32');
PacketSize=fread(fid,1,'uint32');
TimeReslutionTimeStamps=fread(fid,1,'uint32');
TimeReslutionSamples=fread(fid,1,'uint32');
%unsure about actualtype of TimeOrigin
TimeOrigin=fread(fid,8,'uint16');
Application=fread(fid,32,'*char');
Comment=fread(fid,256,'*char');
ExtendedHeaderNumber=fread(fid,1,'ulong');

%Waitbar ------------------------------------------------------------------
waitbar (1/5);
%--------------------------------------------------------------------------

%read extended headers
for i=1:ExtendedHeaderNumber
    Identifier=fread(fid,8,'*char')';
    %modify this later
    switch Identifier
        case 'NEUEVWAV'
            ElecID=fread(fid,1,'uint16');
            PhysConnect=fread(fid,1,'uchar');
            PhysConnectPin=fread(fid,1,'uchar');
            FileInfo.nVperBit(ElecID)=fread(fid,1,'uint16');
            EnergyThresh=fread(fid,1,'uint16');
            FileInfo.HighThresh(ElecID)=fread(fid,1,'int16');
            FileInfo.LowThresh(ElecID)=fread(fid,1,'int16');
            SortedUnits=fread(fid,1,'uchar');
            BytesPerSample=((fread(fid,1,'uchar'))>1)+1;
            temp=fread(fid,10,'uchar');
        otherwise, % added26/7/05 after identifying bug in reading etended headers
            temp=fread(fid,24,'uchar');
    end
end

%Waitbar ------------------------------------------------------------------
waitbar (2/5);
%--------------------------------------------------------------------------

% Calculate number of packets
fseek(fid,0,1);
FileSize=ftell(fid);
PacketNumber=(FileSize-HeaderSize)/PacketSize;

%initialization
fseek(fid,HeaderSize,-1);

WaveformSize=(PacketSize-8)/BytesPerSample;

ByteLength=[num2str(WaveformSize) '*int' num2str(BytesPerSample*8) '=>int' num2str(BytesPerSample*8)];

times = double(fread(fid, PacketNumber, '*uint32', PacketSize-4));
fseek(fid, HeaderSize+4, -1);
chans = uint8(fread(fid, PacketNumber, '*uint16', PacketSize-2));
fseek(fid, HeaderSize+8, -1);
spikes = fread(fid, [WaveformSize PacketNumber], ByteLength, 8)';

fclose(fid);

%Waitbar ------------------------------------------------------------------
waitbar (3/5);
%--------------------------------------------------------------------------

SpikesNumber = zeros(1,255);
for i1=1:max(chans)
    SpikesNumber(i1)=length(find(chans==i1));
end
ActiveChannels=find(SpikesNumber);

ShockThreshold=round(0.9*length(ActiveChannels));
suspects=find((times(ShockThreshold:end)-times(1:end-ShockThreshold+1))<=2);
for i2=suspects'
    ShockPackets= find(abs(times-times(i2))<2);
    chans(ShockPackets)=0;
end;

S = cell(length(ActiveChannels), 1);

%Waitbar ------------------------------------------------------------------
waitbar (4/5);
%--------------------------------------------------------------------------

for i = 1:length(ActiveChannels)
    S{i}.Channel = ActiveChannels(i);
    
    idx = (chans == ActiveChannels(i));
    chans(idx) = [];
    S{i}.Times = times(idx);
    times(idx) = [];
    S{i}.Spikes = spikes(idx, :);
    spikes(idx, :) = [];
end;

%Waitbar ------------------------------------------------------------------
waitbar (5/5);
close (Wb);
%--------------------------------------------------------------------------
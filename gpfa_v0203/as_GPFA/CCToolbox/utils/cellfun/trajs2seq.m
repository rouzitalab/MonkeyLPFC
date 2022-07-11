function [y,x,Seq] = trajs2seq(trajs,zero,flen,DoSeq)
%Trajs2Seq  Convert 'Trajs' struct to the Sequence (or Feature Vector) format.
%
%   See HELP CCToolbox for an explanation of the 'Trajs' argument, the Sequence 
%   format, the Feature Vector format, the Cell format, and all other
%   general concepts referred herein.
%
%   [Y,x,Seq] = Trajs2Seq(TRAJS,[ZERO],[FIXED_LEN])
%   This is the standard form. It takes a set of curves in TRAJS and 
%   outputs the curves in Sequence format. The use of ZERO is given below. FIXED_LEN, if
%   non-zero, restricts the output curves to have the same length
%   given by FIXED_LEN. 
%
%   [Y,DIMS] = Trajs2Seq(TRAJS,[ZERO],[FIXED_LEN],'matrix')
%   This is an alternate form provided for historical reasons. This 
%   performs the same as above except that the curves are returned in 
%   Feature Vector format. If FIXED_LEN is not specified, then all output 
%   curves will have the same length as the shortest input curve. 
%   DIMS gives the number of dimensions of the input curves.
%
%   zero
%     - 'zero' : first measurement is subtracted from each trajectory
%     - 'zero_nocut' : first zeroed entry is left in the vector
%     - 'mean' : trajectory mean is subtracted
%     - 'norm' : trajectory mean is subtracted and then the standard deviation
%                is divided into the trajectory (normalized)
%     - 'znorm': first measurement is subtracted from a normalized trajectory
%     - 'znorm_nocut': first zeroed entry is left in the vector
%     otherwise: if zero is not one of these strings, then the
%                trajectory is left un-transformed. 

% Scott J Gaffney   10 December 2002
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% ---------------------------------
%

PROGNAME = 'trajs2seq';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


NOTHING      = 0;
FIRST        = 1;
MEAN         = 2;
NORM         = 3;
ZNORM        = 4;

Seq=[];   % removes 'output variable unassigned' errors
DoCut = 0;
zero = cexist('zero','');
flen = cexist('flen',0);
DoSeq = cexist('DoSeq',1);

%%% Handle Argument Processing
%%%
% trajs conversion
if (isfield(trajs,'Seq'))
  y = trajs.Y;  x = trajs.X;  Seq = trajs.Seq;
  return;
end
xdesign = [];  % used for fixed design intervals (or non-random x-values)
if (~isfield(trajs,'Y')), trajs.Y = trajs; end
if (~isfield(trajs,'X')), trajs.X = [];    end
if (~iscell(trajs.X))
  xdesign = trajs.X(:)';
  trajs.X = [];
end

% output format
if (isstr(DoSeq) & strcmp(DoSeq,'matrix'))
  DoSeq = 0;
end  

% possible transformations
switch (zero)
  case 'zero'
    zero = FIRST;
    DoCut = 1;
  case 'zero_nocut'
    zero = FIRST;
  case 'mean'
    zero = MEAN;
  case 'norm'
    zero = NORM;
  case 'znorm'
    zero = ZNORM;
    DoCut = 1;
  case 'znorm_nocut'
    zero = ZNORM;
  otherwise
    zero = NOTHING;
end
%%%
%%% End Argument Processing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert if curves are in a cell array
if (iscell(trajs.Y))
  n = length(trajs.Y);
  [mn, lens] = meanlength(trajs.Y);
  minlen = min(lens);
  if (flen>minlen)
    error(['%s: some trajectories don''t meet requested ',...
        'size of %d.\n'], PROGNAME,flen);
  end
  if (isempty(trajs.X))
    if (isempty(xdesign))
      xdesign = (0:max(lens)-1);
    elseif (length(xdesign)<max(lens))
      error(['%s: specified x-design does not meet maximum length', ...
          ' of %d for curves.\n'], PROGNAME,max(lens));
    end
  end
  
  N = sum(lens);
  if (DoSeq & flen~=0), N=n*flen; end
  if (~DoSeq & flen==0), flen=minlen; end
  D = size(trajs.Y{1},2);
  
  if (DoSeq)
    Seq = ones(1,n+1);
    if (DoCut)
      x = zeros(N-n,1);
      y = zeros(N-n,D); 
    else  
      x = zeros(N,1);
      y = zeros(N,D); 
    end
  else
    x = D;      % feature vectors have no independent variables, so pass...
    if (DoCut)  % ...the dimension as the second output variable instead.
      y = zeros(n,D*flen-D);  % we remove the first constant
    else
      y = zeros(n,D*flen);
    end
  end

  % copy the data
  for i=1:n
    tmp = trajs.Y{i};
    if (flen~=0), ni=flen;  else ni=size(tmp,1); end
    % handle transformations
    if (zero==FIRST)
      tmp = tmp(1:ni,:) - ones(ni,1)*tmp(1,:);
    elseif (zero==MEAN)
      tmp = tmp(1:ni,:) - ones(ni,1)*mean(tmp);
    elseif (zero==NORM)
      tmp = (tmp(1:ni,:)-ones(ni,1)*mean(tmp))./(ones(ni,1)*std(tmp));
    elseif (zero==ZNORM)
      tmp = (tmp(1:ni,:)-ones(ni,1)*tmp(1,:))./(ones(ni,1)*std(tmp));
    else  % NOTHING
      tmp = tmp(1:ni,:);
    end
    startx = 1;
    leni = ni;
    if (DoCut)
      tmp = tmp(2:end,:);  % remove constant from beginning
      startx = 2;
      leni = leni-1;
    end
    if (DoSeq)
      Seq(i+1) = Seq(i)+leni;
      if (~isempty(xdesign))
        x(Seq(i):Seq(i+1)-1) = xdesign(startx:ni);
      else
        x(Seq(i):Seq(i+1)-1) = trajs.X{i}(startx:ni);
      end
      y(Seq(i):Seq(i+1)-1,:) = tmp;
    else
      y(i,:) = tmp(:)';  % concatenate dimensions
    end
  end % copy the data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert if curves are in a feature matrix
else
  [n,m,D] = size(trajs.Y);
  if (flen>m)
    error(['%s: some trajectories don''t meet requested ',...
        'size of %d.\n'], PROGNAME,flen);
  end
  if (~flen), flen=m; end  % set flen to maximum extent
  
  % set up cut/no_cut lengths
  if (DoCut)
    clen=flen-1; 
    startx = 2;
    N = n*clen;
  else 
    clen=flen; 
    startx = 1;
    N = n*flen;
  end

  % set up y and (for sequence format only) build x 
  if (DoSeq)
    Seq = cumsum([1 clen*ones(1,n)]);
    if (isempty(trajs.X))
      if (isempty(xdesign))
        x = (startx-1:flen-1)'*ones(1,n);
      else
        xdesign = xdesign(startx:flen);
        x = xdesign(:)*ones(1,n);
      end
    else
      x = zeros(N,1);
      for i=1:n
        x(Seq(i):Seq(i+1)-1) = trajs.X{i}(startx:flen);
      end
    end
    x=x(:);
    y = zeros(N,D);
  else
    x = D; % see comment above for the identical line
    y = zeros(n,clen*D);
  end

  % copy the data
  for d=1:D
    yy = trajs.Y(:,:,d);
    % handle transformations
    if (zero==FIRST)
      yy = yy(:,1:flen) - yy(:,1)*ones(1,flen);
    elseif (zero==MEAN)
      yy = yy(:,1:flen) - mean(yy,2)*ones(1,flen);
    elseif (zero==NORM)
      yy = (yy(:,1:flen)-mean(yy,2)*ones(1,flen))./(std(yy,0,2)*ones(1,flen));
    elseif (zero==ZNORM)
      yy = (yy(:,1:flen)-yy(:,1)*ones(1,flen))./(std(yy,0,2)*ones(1,flen));
    else  % NOTHING
      yy = yy(:,1:flen);
    end
    if (DoCut),  yy = yy(:,2:end); end
    if (DoSeq)
      yy = yy';  
      y(:,d) = yy(:);
    else
      indx = (d-1)*clen+1:d*clen;
      y(:,indx) = yy;
    end
  end  
end  % handle feature vector
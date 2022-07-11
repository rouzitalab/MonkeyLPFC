function [samp, state] = sampfun(fun, N, stdv, range, seed)
%SAMPFUN  Draw a sample of points from a polynomial function.
%   SAMP = SAMPFUN(FUN,N) is an N by 2 matrix containing x,y pairs
%   of coordinates for each of the N points that are "drawn" from the
%   function FUN.  FUN is a vector similar to that taken by ROOTS.
%
%   SAMP = SAMPFUN(FUN,N,STDV) is an N by 2 matrix like that above, where
%   the y-coordinates for each corresponding x-coordinate are
%   tainted with Gaussian noise of mean 0 and standard deviation STDV.
%
%   SAMP = SAMPFUN(FUN,N,STDV,RANGE) is an N by 2 matrix like that above,
%   where the N sample points are drawn uniformly from the range
%   RANGE.
%
%   [SAMP, STATE] = SAMPFUN(FUN,N,STDV,RANGE,SEED) returns an N by 2
%   matrix like that above as SAMP, and returns the current state of
%   the uniform psuedo-random number generator as STATE.
%   You may set the seed (or the state) of the uniform random number
%   generator by giving SEED (e.g., use STATE from the previous call). 
%
%   NOTE: if RANGE is not specified (or is empty []), then the range 
%   will be set to the "smallest" interval that still includes all of 
%   the real roots of FUN as determined by ROOTS(FUN). If SEED is
%   empty or is not specified, then the current state will be used.

% Scott Gaffney   26 April 1998
% Department of Information and Computer Science
% University of California, Irvine

% Handle parameter passing processing
if (nargin < 2)
    error('error: you must supply at least 2 arguments; see help sampfun');
end
if (isempty(fun) | ndims(fun) > 2 | min(size(fun)) ~= 1)
    error('error: FUN must be a nonempty vector');
end
fun = shiftdim(fun);  % make a col vector
if (isempty(N) | prod(size(N)) ~= 1)
    error('error: N must be a scalar');
end
if (exist('stdv', 'var') ~= 1)
    stdv = 0;
end
if (exist('range', 'var') ~= 1 | isempty(range))
    realrange = [];
    rts = roots(fun);
    for i = 1:length(rts)
	if (isreal(rts(i)))
	    realrange = [realrange; rts(i)];
	end
    end
    if (isempty(rts) | isempty(realrange) | length(realrange) < 2)
	error('error: FUN does not have at least two real roots');
    end
    range = [(min(realrange)-100*eps) (max(realrange)+100*eps)];
end

% Select x-coords uniformly and evaluate FUN
if (exist('seed', 'var') == 1 & ~isempty(seed))
    rand('state', seed);
end
state = rand('state');
samp = zeros(N,2);
samp(:,1) = rand(N,1)*(range(2)-range(1)) + range(1);
samp(:,2) = polyval(fun, samp(:,1)) + randn(N,1)*stdv;


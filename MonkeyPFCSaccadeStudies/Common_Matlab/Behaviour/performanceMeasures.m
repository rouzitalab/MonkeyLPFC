function [dprime, bias, accuracy] = performanceMeasures(TrueA, FalseA, FalseB, TrueB, varargin)

%If probability approaches 0 or 1, fix by pmarg
if nargin > 4
    pmarg = varargin{1};
else
    pmarg = 0.01;
end

% Calculate performance measures
%         cnfMat = [TrueA FalseA; ...
%                   FalseB TrueB];

%http://www.birmingham.ac.uk/Documents/college-les/psych/vision-laboratory/sdtintro.pdf
tpr = TrueA./(TrueA + FalseB);
tpr(tpr<pmarg)=pmarg; tpr(tpr>1-pmarg)=1-pmarg;
fpr = FalseA./(FalseA+TrueB);
fpr(fpr<pmarg)=pmarg; fpr(fpr>1-pmarg)=1-pmarg;
ztpr = norminv(tpr, 0, 1);
zfpr = norminv(fpr, 0, 1);
dprime = ztpr - zfpr;
bias = -(ztpr+zfpr)/2;

%70% for both classes is d=1.0488.
%Highest possible is 6.93, but effectively 4.65 for 99%

%Other measures of performance.
% sens = tp ./ (tp+fp);
% spec = tn ./ (tn+fn);
% balAcc = (sens+spec)/2;
% informedness = sens+spec-1;

accuracy = 100* (TrueA + TrueB) / (TrueA + FalseA + FalseB + TrueB);
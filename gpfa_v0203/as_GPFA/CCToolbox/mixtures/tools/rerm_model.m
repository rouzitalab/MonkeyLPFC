function Model = rerm_model(B,Sigma,Alpha,Gamma,R,Pik,Lhood,C,NumPoints,...
    back,order,zero,method)
%RERM_MODEL  Makes a model structure for random effects regression mixtures
%   Model = RERM_MODEL(B,Sigma,Alpha,Gamma,R,Pik,Lhood,C,NumPoints, ...
%       back,order,zero,method)
%

% Scott Gaffney   05 April 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'rerm_model';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


Model.B = B;
Model.Sigma = Sigma;
Model.Alpha = Alpha;
Model.Gamma = Gamma;
Model.R = R;
Model.Pik = Pik;
Model.Lhood = Lhood;
Model.C = C;
Model.NumPoints = NumPoints;
Model.back = cexist('back',0);
Model.order = cexist('order',[]);
Model.zero = cexist('zero','nozero');
Model.method = cexist('method','rerm');

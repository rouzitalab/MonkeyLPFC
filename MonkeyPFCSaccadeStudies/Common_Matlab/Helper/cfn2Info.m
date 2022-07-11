function [pCorr, infPerTrial, H] = cfn2Info(cfnMat)
%[pCorr, infPerTrial, H] = cfn2Info(cfnMat)
nClasses = size(cfnMat,1);
pCorr = sum(diag(cfnMat)) / sum(sum(cfnMat));
p0 = 1/nClasses;
infPerTrial = sign(pCorr-p0) * ( pCorr*log2(pCorr/p0) + (1-pCorr)*log2( (1-pCorr)/(1-p0) ) );
Ntot = sum(sum(cfnMat));
cfnMat(cfnMat==0)=10^(-10);
A=0;
[I,J] = size(cfnMat);
for i1 = 1:I
    for j1 = 1:J
        B=0;
        C=0;
        for i2 = 1:I
            B=B+cfnMat(i2,j1);
        end
        for j2 = 1:J
            C=C+cfnMat(i1,j2);
        end
        A=A+ cfnMat(i1,j1)*[log2(cfnMat(i1,j1)) - log2(C) - log2(B) + log2(Ntot)];
    end
end
H= (1/Ntot)*A;
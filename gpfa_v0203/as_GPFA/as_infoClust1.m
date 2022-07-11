% attempt to recreate Victor and Purpura info of clustering

function H = as_infoClust1(N);

% N is confusion matrix, and H is info
% i goes to I
% j goes to J

Ntot = sum(sum(N));
 N(find(N==0))=.0000000001;
A=0;
[I J] = size(N);
    
for i1 = 1:I
    for j1 = 1:J
        B=0;
        C=0;
        for i2 = 1:I
            B=B+N(i2,j1);
        end
        for j2 = 1:J
            C=C+N(i1,j2);
        end
        A=A+ N(i1,j1)*[log2(N(i1,j1)) - log2(C) - log2(B) + log2(Ntot)];
    end
end

H= (1/Ntot )*A;


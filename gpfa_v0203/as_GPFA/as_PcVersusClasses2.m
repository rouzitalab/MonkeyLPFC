% needs further work since real classes you can cheat and add adjacents,
% but resorted classes need to have flexibility (ie use factorial rule)
% this assumes that each grouping of n classes is the superset for each
% grouping of n-1 classes
% written ajs Aug 2011

s=size(cm,2);
w=[];
PC=[];
 cm1= cm;
for i = s:-1:3  % ie for every possibly number of categories
    t=size(cm1,2);
    Pcarray=[];
    canIdx=0;
    canIdxLookup=[];
    for j = 1:t % figure out which row to eliminate
        can_cm2 = cm1;
        if j==t
            can_cm2(1,:)= can_cm2(1,:)+can_cm2(j,:);
            can_cm2(j,:)=[];
        else
            can_cm2(j+1,:)= can_cm2(j+1,:)+can_cm2(j,:);
            can_cm2(j,:)=[];
        end
        for k = 1:t % figure out which column to eliminate
            for l=1:t % figure out which column to add it to
                if l~=k
                    can_cm3=can_cm2;
                    canIdx=canIdx+1;
                    canIdxLookup = [canIdxLookup;canIdx j k l];
                    can_cm3(:,l)=can_cm3(:,l)+can_cm3(:,k);
                    can_cm3(:,k)=[];
                    Pcarray(canIdx) = sum(diag(can_cm3)) / sum(sum(can_cm3));
                end
            end
        end
    end
    [PcVal jmax] = max(Pcarray);
    PC(i-2)=PcVal;
    [i PcVal];
    J=canIdxLookup(jmax,2);
    K=canIdxLookup(jmax,3);
    L=canIdxLookup(jmax,4);
    w(i).cm2 = cm1;
    if J==size(cm1,2)
        w(i).cm2(1,:) = cm1(1,:)+cm1(J,:);
    else
    w(i).cm2(J+1,:) = cm1(J+1,:)+cm1(J,:);
    end
    w(i).cm2(:,L) = cm1(:,L)+cm1(:,K);
  w(i).cm2(J,:)=[];
   w(i).cm2(:,K)=[];
   cm1=w(i).cm2;
    
end
PC=[Pc fliplr(PC)];
 figure,plot([8:-1:2],PC),axis([2 8 min(PC)-.1 1.1]),set(gca,'xdir','rev')
        
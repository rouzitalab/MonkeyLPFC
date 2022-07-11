% needs further work since real classes you can cheat and add adjacents,
% but resorted classes need to have flexibility (ie use factorial rule)

s=size(cm,2);

PC=[];
 cm1= cm;
for i = s:-1:2
    t=size(cm1,2);
    Pcarray=[];
    for j = 1:t-1
        can_cm2 = cm1;
        can_cm2(j+1,:)= can_cm2(j+1,:)+can_cm2(j,:);
        can_cm2(:,j+1)=can_cm2(:,j+1)+can_cm2(:,j);
        can_cm2(j,:)=[];
        can_cm2(:,j)=[];
        
        Pcarray(j) = sum(diag(can_cm2)) / sum(sum(can_cm2));
    end
    [PcVal jmax] = max(Pcarray);
    PC(i)=PcVal;
    k(i).cm2 = cm1;
    k(i).cm2(jmax+1,:) = cm1(jmax+1,:)+cm1(jmax,:);
    k(i).cm2(:,jmax+1) = cm1(:,jmax+1)+cm1(:,jmax);
   k(i).cm2(jmax,:)=[];
   k(i).cm2(:,jmax)=[];
   cm1=k(i).cm2;
    
end
PC=fliplr(PC);
        
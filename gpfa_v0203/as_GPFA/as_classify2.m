function [orderedGroupEst,map, confmat] = as_classify2(origidx,unorderedGroupEst)

% requires groups to be named from 1 to max
% uses origidx, instead of origGroup
% eg origGroup{1}=[1 1 1 1 1 1];
n=max(origidx);
for i=1:n
    origGroup{i}=i*ones(sum(origidx==i),1);
    origGroupIndex{i} = find(origidx==i);
end
uGE = unorderedGroupEst;
G = origGroup;
% oGE = orderedGroupEst;
maps = perms(1:numel(G));
k = size(maps,1);
dist=zeros(k,1);
SUMS=zeros(k,1);
tic
for i=1:k
    Gmeans=[];
    uGEmeans=[];
    sumofmatches=[];
    for j=1:numel(G)
        Gmeans(1,j) = mean(G{j});
        uGEmeans(1,j) = mean(uGE{maps(i,j)});
        a=uGE{maps(i,j)};
        b=G{j};
        sumofmatches(1,j)= sum(ismember(a,b));
    end
    SUMS(i) = sum(sumofmatches);
    dist(i)=as_Euclid(Gmeans,uGEmeans);
   % clear sumofmatches
end
 toc       
%[mini mapidx]=min(dist);
[maxi mapidx]=max(SUMS);
map=maps(mapidx,:);

confmat=zeros(numel(G));

for i=1:numel(G)
     oGE{i} = uGE{map(i)};
     for j=1:numel(G)
         confmat(i,j) = sum(oGE{i}==j);
     end
end

orderedGroupEst = oGE;



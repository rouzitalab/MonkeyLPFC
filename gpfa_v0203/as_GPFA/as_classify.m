function [orderedGroupEst,map, confmat] = as_classify(origGroup,unorderedGroupEst)

% eg origGroup{1}=[1 1 1 1 1 1];
dumbbell=0;

uGE = unorderedGroupEst;
% no see below - G = origGroup;
% step 1 get rid of aPeriods with no trials
S=size(origGroup,2);
S2=size(uGE,2);
if dumbbell
    for s=1:S/2
        origGroup{s}=s*ones(size(origGroup{s},1)+size(origGroup{s+4},1),1)
        origGroup{s+4}=[];
    end
    for s=1:S2
        ind=uGE{s}>4;
        uGE{s}(ind)=uGE{s}(ind)-4;
    end
        
end
        
g=1;
for s=1:S
    if ~isempty(origGroup{s})
        G{g}=g*ones(size(origGroup{s}));
        g=g+1;
    end
end
% oGE = orderedGroupEst;
maps = perms(1:numel(G));
k = size(maps,1);
dist=zeros(k,1);
SUMS=zeros(k,1);
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
    clear sumofmatches
end
        
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



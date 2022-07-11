function [orderedGroupEst,map, confmat] = as_classify3(origGroup,unorderedGroupEst,classMethod)

% eg origGroup{1}=[1 1 1 1 1 1];
%classify 3 uses max number of matching numbers (instead of distance
%between means). May also expand to allow option to optimize Pc and H
% classMethod =1 : match total number of correct classifications in all
% groups
% classMethod =2 : optimize Pc
% classMethod =3 : optimize H

if nargin ==2
    classMethod = 1;
end

dumbbell=0;


uGE = unorderedGroupEst;
% no see below - G = origGroup;
% step 1 get rid of aPeriods with no trials
S=size(origGroup,2);
S2=size(uGE,2);
if dumbbell
    for s=1:S/2
        origGroup{s}=s*ones(size(origGroup{s},1)+size(origGroup{s+4},1),1);
        origGroup{s+4}=[];
    end
    for s=1:S2
        ind=uGE{s}>4;
        uGE{s}(ind)=uGE{s}(ind)-4;
    end
        
end
   
alphabet=['abcdefghijklmnopqrstuvwxyz'];
lettersUsed=[];
for s=1:S2
    N=size(uGE{s},2);
    for n=1:N
        
        val = uGE{s}(n);
        uGE2{s}(n)=alphabet(val);
        lettersUsed=[lettersUsed alphabet(val)];
    end
end
lettersUsed=unique(sort(lettersUsed));



realGroupNames=[];
g=1;
for s=1:S
    if ~isempty(origGroup{s})
        G{g}=g*ones(size(origGroup{s}));
        realGroupNames=[realGroupNames g];
        eval([lettersUsed(g) '= g' ]); 
        g=g+1;
    end
end

%now reassign the group names so they are comparable
for s=1:S2
    N=size(uGE{s},2);
    for n=1:N
        uGE3{s}(n)=eval(uGE2{s}(n));
    end
end

% clear uGE
uGE=uGE3;

% oGE = orderedGroupEst;
maps = perms(realGroupNames);
% ie maps(1,:)=[3 2 1]
% this should be read as "the third real group is mapped to the 
% first cluster group
NumClusters = numel(uGE);
if NumClusters < size(maps,2)  % in case there are less clusters than original groups
    maps = maps(:,NumClusters);
end
k = size(maps,1);
%dist=zeros(k,1);
%SUMS=zeros(k,1);
CategoriesAssigned=zeros(1,numel(G));
totalMatches=[];
for i=1:k
    Gmeans=[];
    uGEmeans=[];
    sumofmatches=[];
    for j=1:numel(G)
        % Gmeans(1,j) = mean(G{j});
        % uGEmeans(1,j) = mean(uGE{maps(i,j)});
        a=uGE{maps(i,j)};
        b=G{j};
        sumofmatches(1,j)= sum(ismember(a,b));
        if classMethod>1
            for ii=1:numel(G)
                oGE{ii} = uGE{maps(i,ii)};
                for j=1:numel(G)
                    confmat(ii,j) = sum(oGE{ii}==j);
                end
            end
            % clear oGE
        end
    end
    if classMethod>1
        PcArray(1,i) = sum(diag(confmat)) / sum(sum(confmat));
        HArray(1,i) = as_infoClust1(confmat);
    end
    SUMS(i) = sum(sumofmatches);
    dist(i)=as_Euclid(Gmeans,uGEmeans);
    % clear sumofmatches
end

%[mini mapidx]=min(dist);
switch classMethod
    case 1
        [maxi mapidx]=max(SUMS);
    case 2
        [maxi mapidx]=max(PcArray);
    case 3
        [maxi mapidx]=max(HArray);
end
[maxi mapidx]=max(SUMS);
map=maps(mapidx,:);
CategoriesAssigned(1,map)=1;
j=1;

while sum(CategoriesAssigned)<length(CategoriesAssigned)  %NOT YET WORKING
    
    ZeroInds=find(CategoriesAssigned==0);
    maps=perms(ZeroInds);
    
    FirstInd=ZeroInds(1);
    a=uGE{maps(i,j)};
    b=G{j};
    sumofmatches(1,j)= sum(ismember(a,b));
    if NumClusters < size(maps,2)  % in case there are less clusters than original groups
        maps = maps(:,NumClusters);
    end
    k = size(maps,1);
    
end

confmat=zeros(numel(G));

for i=1:numel(G)
    oGE{i} = uGE{map(i)};
    for j=1:numel(G)
        confmat(i,j) = sum(oGE{i}==j);
    end
end

orderedGroupEst = oGE;



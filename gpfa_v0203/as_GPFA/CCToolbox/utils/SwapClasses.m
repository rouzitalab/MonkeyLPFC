function model = SwapClasses(model,truemodel)
%SwapClasses  Swap classes/fields according to map or true model.
%   new_model = SwapClasses(model,truemodel)
%
%   new_model = SwapClasses(map,model) uses the map as
%   follows: new_model.C(find(model.C==map(i))) = i.

% Scott Gaffney   08 June 2002
% Department of Information and Computer Science
% University of California, Irvine


PROGNAME = 'SwapClasses';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

% mapsig = cexist('mapsig','mapsigma');
% if (strcmp(mapsig,'mapsigma'))
%   mapsig = 1;
% else
%   mapsig = 0;
% end

if (~isstruct(model))
  map = model;
  model = truemodel;
  sortc = zeros(length(model.C),1);
  uniq = unique(model.C);
  for i=uniq(:)'
    sortc(find(model.C==i)) = find(map==i);
  end
else
  [sortc,map] = sortcls(model.C,truemodel.C);
end


%% Alpha
if (isfield(model,'Alpha'))
  model.Alpha = model.Alpha(map);
end


%% Mu
if (isfield(model,'Mu'))
  model.Mu = model.Mu(:,map,:);
end


%% Sigma
if (isfield(model,'Sigma'))
  if (ndims(model.Sigma)==2)
    if (isvector(model.Sigma))
      model.Sigma = model.Sigma(map);
    else
      model.Sigma = model.Sigma(map,:);
    end
  else
    model.Sigma = model.Sigma(:,:,map,:);
  end
end


%% R
if (isfield(model,'R'))
  if (strcmp(model.method,'rerm'))  % not really supported anymore
    if (ndims(model.R)==2)
      model.R = model.R(map,:);
    else
      model.R = model.R(:,:,map,:);
    end
  else
    if (isvector(model.R))
      model.R = model.R(map);
    else
      model.R = model.R(map,:);
    end
  end
end


%% S
if (isfield(model,'S'))
  if (isvector(model.S))
    model.S = model.S(map);
  else
    model.S = model.S(map,:);
  end
end


%% U
if (isfield(model,'U'))
  if (isvector(model.U))
    model.U = model.U(map);
  else
    model.U = model.U(map,:);
  end
end


%% T
if (isfield(model,'T'))
  if (isvector(model.T))
    model.T = model.T(map);
  else
    model.T = model.T(map,:);
  end
end


%% Pik
if (isfield(model,'Pik'))
  model.Pik = model.Pik(:,map);
end


%% Ea
if (isfield(model,'Ea'))
  model.Ea = model.Ea(:,map,:);
  model.Va = model.Va(:,map,:);
  if (isfield(model,'Vab'))
    model.Vab = model.Vab(:,map,:);
  end
  if (isfield(model,'Vae'))
    model.Vae = model.Vae(:,map,:);
  end
  if (isfield(model,'Vaf'))
    model.Vaf = model.Vaf(:,map,:);
  end
end


%% Eb
if (isfield(model,'Eb'))
  model.Eb = model.Eb(:,map,:);
  model.Vb = model.Vb(:,map,:);
  if (isfield(model,'Vbe'))
    model.Vbe = model.Vbe(:,map,:);
  end
  if (isfield(model,'Vbf'))
    model.Vbf = model.Vbf(:,map,:);
  end
end


%% Ee
if (isfield(model,'Ee'))
  model.Ee = model.Ee(:,map,:);
  model.Ve = model.Ve(:,map,:);
  if (isfield(model,'Vef'))
    model.Vef = model.Vef(:,map,:);
  end
end


%% Ef
if (isfield(model,'Ef'))
  model.Ef = model.Ef(:,map,:);
  model.Vf = model.Vf(:,map,:);
end




model.map = map;
model.unsortC = model.C;
model.C = sortc;

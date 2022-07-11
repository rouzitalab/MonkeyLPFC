function H = pick(i)
%PICK  Toggles the selection of a curves in the current plot
%  H = PICK(i) toggles the selection of the curve(s) with Tag property set to
%  the string num2str(i) and returns its handle.

% Scott J. Gaffney  5 February 2004
% School of Computer Science
% University of California, Irvine


PROGNAME = 'pick';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


h=[];
for j=1:length(i)
  f = findobj(gca,'Tag',num2str(i(j)));
  if (isempty(f)), continue; end
  h(end+1) = f;
  select = get(h(end),'selected');
  if (strcmp(select,'on'))
    set(h(end),'selected','off');
  else
    set(h(end),'selected','on');
  end
end
if (nargout==1)
  H=h;
end
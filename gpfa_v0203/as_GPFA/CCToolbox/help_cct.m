function help_cct(x)
%Help_CCT   Starts the CCToolbox help system using the matlab helpwin browser.
%   Help_CCT() will open the matlab helpwin browser and display the
%   root Contents.m file of the CCToolbox.
%
%   Help_CCT(ShowInCommandWindow) will show the root Contents.m file
%   in the matlab command window instead of in the helpwin browser.
%   ShowInCommandWindow can be any value as long as it exists (e.g., it can
%   even be empty).

%  Scott J. Gaffney   18 May 2005
%  University of California, Irvine


base = fileparts(mfilename('fullpath'));

if (exist('x')~=1)
  helpwin(base);
else
  more on;
  help(base)
  more off;
end
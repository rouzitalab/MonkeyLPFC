% Cell Array Functions (CCToolbox/cellfun)
%  CCToolbox 1.0
%
% Data format conversion
%  trajs2seq    - Convert 'Trajs' struct curves to Sequence format
%  trajs2cell   - Convert 'Trajs' struct curves to Cell format
%  seq2cell     - Convert Sequence format curves to Cell Format
%  fvec2cell    - Convert Feature-vector curves to Cell format
%  trajslen     - Find number of curves in 'Trajs' struct
%  cell2nan     - Convert cell array into a matrix with NaNs for missing values
%  cell2vector  - Flatten-out a cell array
%
% Cell Statistics
%  cellmax      - Max of each cell
%  cellmin      - Min of each cell
%  cellmedian   - Median of each cell
%  cellmean     - Mean of each cell
%  cellstd      - Standard deviation of each cell
%  cell_len     - Length of each cell
%  meanlength   - Mean of all cell lengths, equal to mean(cellmean(X))
%
% Other
%  trajs2regparam - (Out-dated) Convert curves to regression parameters
%
% See also CCTOOLBOX

% Scott J. Gaffney   18 May 2005
% University of California, Irvine
% The Curve Clustering Toolbox 1.0
%  Scott J. Gaffney   18 May 2005
%  University of California, Irvine
%
% The Curve Clustering Toolbox is a Matlab toolbox that implements a family
% of probabilistic model-based curve-aligned clustering algorithms. The 
% cluster models themselves are based on polynomial and spline regression 
% mixture models that allow for continuous curve alignment in both
% measurement space and in time. Learning is carried out using an
% EM (Expectation-Maximization) framework. The model specification and
% learning framework are detailed in (Gaffney, 04) which can be found 
% online at http://www.ics.uci.edu/~sgaffney/papers/sgaffney_thesis.pdf
%
% Important: You must run SetCCTPath() each time you wish to use this
% toolbox. You can run this automatically by starting matlab with the
% -r option.
%
% This toolbox is organized within various directories as follows:
%  CCToolbox\docs\             - Documentation
%  CCToolbox\gui\              - Graphical user interface tools
%  CCToolbox\mixtures\         - Clustering models
%  CCToolbox\model_selection\  - Model selection tools
%  CCToolbox\regression\       - Regression tools
%  CCToolbox\simulation\       - Data simulation tools
%  CCToolbox\splines\          - Spline and spline regression tools
%  CCToolbox\stats\            - Some statistics functions
%  CCToolbox\utils\            - General utilities
%  CCToolbox\view\             - Visualization tools
%
% You can view this help file from within the matlab helpwin browser by
% running help_cct() at any time. You can get a list of the functions in each 
% of the above directories by typing 'helpwin <directory>' or 'help <directory>' 
% at the matlab prompt. If you are currently viewing this help file from within 
% matlab's helpwin browser, then you can simply click the hyperlinks above as 
% a shortcut.
%
% Much of the clustering and visualization functionality of the CCToolbox can
% be accessed through two simple functions: 
%
%  Curve_Clust  - Cluster trajectories using chosen method
%  ShowModel    - Plot clustering model parameters and/or data
%
% Curve_Clust() takes a set of curves in various supported data structures and 
% runs the desired clustering method, returning the results in a standard model 
% structure. The results of the clustering can be plotted by passing the returned 
% model structure to ShowModel(). See also CCT_CURVE_FORMAT and CCT_MODEL_STRUCT
% for details on the supported data and model structures.
%
% The quickest way to begin using this toolbox is to read CCT_CURVE_FORMAT and
% CCT_RUN_CLUSTERING, and then refer to Curve_Clust() to begin exploratory
% clustering of your curve data. In addition, there are a number of useful 
% documents provided with this toolbox that should read. These 
% documents are located in the docs\ folder. You can get a list of these 
% documents by clicking the below link or by typing 'help cctoolbox\docs' at 
% the matlab prompt.
%
%  CCToolbox\docs\             - Documentation

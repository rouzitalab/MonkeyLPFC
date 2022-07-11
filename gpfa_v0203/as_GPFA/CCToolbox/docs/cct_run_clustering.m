function cct_run_clustering(x)
%CCT_Run_Clustering
% Curve Clustering Toolbox 1.0
%  How to begin clustering
%
%  You can cluster curves by following four simple steps.
%
%  1. Load the appropriate curve data and place the curves into the appropriate
%     'Trajs' format; see CCT_CURVE_FORMAT.
%  2. Set the cluster model options (e.g., number of clusters).
%  3. Set the EM algorithm options (e.g., number of EM starts).
%  4. type 'model = curve_clust(trajs,options)' at the matlab prompt.
%
%  Important note: you can easily determine what default options are used
%  by any cluster model by simply calling the clustering algorithm with
%  the single string argument 'options'. For example, srm_d('options') will 
%  print out the default options used for the SRM_D cluster model.
%
%  Setting cluster model options
%  -----------------------------
%  There are a number of cluster model options. Each model may require
%  different options; e.g., knots are required for spline regression but
%  not for polynomial regression. Listed below are some common required options.
%
%  ops.method         : clustering method. See also LISTMODELS for values
%  ops.K              : number of clusters
%  ops.order          : specify (1) linear, (2) quadratic, etc.
%  ops.zero           : zero data first? See trajs2seq() for values
%  ops.knots          : set of knots for spline regression
%
%  Other options are truly 'optional'. Some of these options are shown below.
%
%  ops.ShowGraphics   : (0,1) enable iterative graphics output during EM
%  ops.Sigma.Share    : (0,1) enable sharing of sigma across D (output dimension)
%  ops.Sigma.Diagonal : (0,1) enable diagonal sigma matrices
%  ops.MinLen         : chop all trajectories to specified length
%
%  Usually, the help message at the beginning of each clustering algorithm
%  will list any unique required options for that model.
%
%  Setting EM algorithm options
%  ----------------------------
%  Besides model options, there are also EM specific options. These are
%  options that are more technical in nature and usually apply to all EM
%  algorithms in general. Some are listed below.
%
%  ops.NumEMStarts    : number of EM starts
%  ops.IterLimit      : maximum number of EM iterations for each start
%  ops.ValidateLhood  : (0,1) recalculate likelihood of the best EM run (this
%                     : can be used when numerical integration is used 
%                     : during EM to get a more accurate estimate of the true
%                     : loglikelihood after all EM starts have been completed).
%  ops.stopval        : value used to terminate EM iterations (see code for info)
%
%  Example curve clustering run
%  ----------------------------
%  As an example, the following four lines will run quadratic polynomial 
%  regression mixtures for the curves in Trajs and find three clusters.
%  Only one EM start will be carried out since one is the default value
%  for ops.NumEMStarts.
%  
%  ops.method = 'lrm';
%  ops.order  = 2;
%  ops.K      = 3;
%  model = curve_clust(Trajs,ops);
%
%  If there are not enough data or the models are too complex, there may be
%  a lot of instances of singular matrices or other warning conditions. 
%  You can stop the display of these warning messages by issuing the
%  following commands.
%
%  warning off all;
%  warning on MATLAB:divideByZero;
%
%  Turning on the divide by zero warning can be useful to detect more
%  serious problems.


%  Scott J. Gaffney   18 May 2005
%  University of California, Irvine

PROGNAME = 'cct_run_clustering';
if (exist('x')~=1)
  helpwin(PROGNAME);
else
  more on;
  help(PROGNAME)
  more off;
end

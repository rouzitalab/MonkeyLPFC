function cct_cluster_models(x)
% CCT_Cluster_Models
%  The Curve Clustering Toolbox 1.0
%   Base Set of Supported Cluster Models
%
%   The current set of supported cluster models can always be listed by running
%   the function ListModels() at the matlab prompt. For example, running this
%   command reveals that K-means is a supported cluster model in this toolbox.
%   You can access the help file for the supported cluster models by following 
%   the CCToolbox\mixtures link from the main CCToolbox help screen.
%   See also CCTOOLBOX.
%
%   There is a specific naming convention that is followed in this toolbox 
%   for the set of regression mixture models. The naming convention is at 
%   first based on the standard linear and spline regression mixture models. 
%   These two models are listed by ListModels() as follows.
%
%      lrm    - Linear Regression Mixtures (LRM) with y = XB + e
%      srm    - Spline Regression Mixture (SRM) with y = XB + e
%
%   This listing shows that the method name 'lrm' should be used to run 
%   LRM and the method name 'srm' should be used for SRM. All of the other 
%   regression mixture models are named by adding a suffix containing the 
%   transformation variables that are learned for that specific model. For 
%   example, the following curve-aligned LRM models are listed by ListModels().
%
%      lrm_d      - LRM with transformation  [x]B + d
%      lrm_b      - LRM with transformation  [x+b]B
%      lrm_d_b    - LRM with transformation  [x+b]B + d
%      lrm_cd     - LRM with transformation c[x]B + d
%      lrm_ab     - LRM with transformation  [ax+b]B
%      lrm_d_ab   - LRM with transformation  [ax+b]B + d
%      lrm_cd_b   - LRM with transformation c[x+b]B + d
%      lrm_cd_ab  - LRM with transformation c[ax+b]B + d
%
%   Notice how the measurement space transformation variables are always listed 
%   before the time transformation variables for each method name (on the 
%   left). A similar but shorter list is returned for the curve-aligned 
%   SRM models.
%
%      srm_d      - SRM with transformation y=[x]B+d
%      srm_cd     - SRM with transformation y=c[x]B+d
%      srm_b      - SRM with transformation y=[x+b]B
%      srm_ab     - SRM with transformation y=[ax+b]B
%      srm_d_b    - SRM with transformation y=[x+b]B+d
%
%   Gaussian mixtures ('gmix') and K-means ('kmeans') are also provided in the
%   base set.

%  Scott J. Gaffney   18 May 2005
%  University of California, Irvine

PROGNAME = 'cct_cluster_models';
if (exist('x')~=1)
  helpwin(PROGNAME);
else
  more on;
  help(PROGNAME)
  more off;
end

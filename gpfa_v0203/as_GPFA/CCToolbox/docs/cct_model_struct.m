function cct_model_struct(x)
%CCT_Model_Struct
% The Curve Clustering Toolbox 1.0
%  Model Structure Description
%
%  Size variables used below for structure definitions
%  ---------------------------------------------------
%   - n   : number of trajectories
%   - D   : number of dimensions (e.g., 2 for lat and lon)
%   - K   : number of clusters
%   - P   : 1 + order of regression model
%   - I   : number of EM iterations
%
%
%  The Base Fields
%  ---------------
%   .K              : (number) number of clusters
%   .order          : (number) order of polynomial regression model
%   .Options        : (struct) stores various option settings
%   .method         : (string) method name, see LISTMODELS
%   .zero           : (string) type of preprocessing, see TRAJS2SEQ
%   .Mu             : (P-K-D matrix) regression parameters
%   .Sigma          : (K-D matrix) variance param for regression noise model
%   .Lhood          : (I-vector) sequence of log-likelihood iterates
%   .NumPoints      : (number) total number of points in all trajs
%   .C              : (n-vector) most likely classification given .Pik
%   .TrainLhood     : (number) equal to .Lhood(end)
%   .TrainLhood_ppt : (number) equal to .TrainLhood/.NumPoints
%   .NumIndParams   : (number) number of independent parameters estimated
%   .state          : (vector) stores RAND state information
%   .nstate         : (vector) stores RANDN state information
%
%
%  The Alignment Fields
%  --------------------
%   .Ea             : (n-K matrix) scale parameters (in time)
%   .Eb             : (n-K matrix) translation parameters (in time)
%   .Ee             : (n-K-D matrix) scale parameters (measurement space)
%   .Ef             : (n-K-D matrix) translation parameters (measurement space)
%   .Va             : (n-K matrix) variance for .Ea
%   .Vb             : (n-K matrix) variance for .Eb
%   .Ve             : (n-K-D matrix) variance for .Ee
%   .Vf             : (n-K-D matrix) variance for .Ef
%   .Vab            : (n-K matrix) covariance for .Ea and .Eb
%   .Vef            : (n-K-D matrix) covariance for .Ee and .Ef
%   .Vae,.Vaf,etc.  : (n-K-D matrix) covariances
%   .Pik            : (n-K matrix) membership probabilities
%   .Alpha          : (K-vector) k-th mixture weight
%   .R              : (K-1 matrix) variance param for prior on .Ea (also for
%                   :  .Ee when .Ea doesn't exist; this should be changed)
%   .S              : (K-1 matrix) variance param for prior on .Eb (also for
%                      .Ef when .Eb doesn't exist; this should be changed)
%   .T              : (K-D matrix) variance param for prior on .Ef
%   .U              : (K-D matrix) variance param for prior on .Ee
%
%
%  Parameter explanations
%  ----------------------
%  In this section, a more detailed description for some of the more complex
%  model structure fields is given.
%
%  As an aside, it is important to note that each of following variables 
%  has an expected value under each cluster (i.e., it's associated structure 
%  field has at least K entries); this is the nature of a mixture model. 
%  However, it is also true that the .E- and .V- parameters are also 
%  trajectory-dependent (i.e., they have at least n*K entries).
% 
%  For example, .Ee gives the expected scale parameter for each trajectory 
%  and each dimension, under each cluster (thus, there are n*D*K entries). 
%  However, the only trajectory-dependent parameter value that is relevant 
%  is the one associated with the corresponding most likely cluster for 
%  each trajectory. We are not normally interested in the expected 
%  transformation for a specific trajectory in a cluster for which 
%  the trajectory does not belong. In other words, we would like a reduced 
%  matrix of size n-by-D that only lists the entries for each trajectory 
%  (and dimension) from the clusters to which the trajectories belong.
% 
%  You can retrieve a reduced matrix that only contains the relevant
%  entries for each trajectory-dependent parameter by using GETCINDX.
%  For example, GETCINDX(model.Ee, model.C) will return an n-D matrix 
%  containing all the relevant entries for all n trajectories and D
%  dimensions (this is explained more clearly for each case below).
%  You can substitute any of the trajectory-dependent 
%  parameters (.Ee, .Ef, .Ve, .Vf, .Vef) for .Ee in the call to
%  GETCINDX as shown below.
%
%  - .Mu   : This variable gives the cluster-dependent regression
%   coefficients. .Mu(:,k,d) gives the coefficients for the k-th cluster
%   in the d-th dimension; note that .Mu(1,k,d) gives the y-intercept.
%
%  - .Sigma : This variable gives the variance parameter of the regression
%   noise model. .Sigma(k,d) gives the parameter for the k-th cluster
%   in the d-th dimension.
%
%  - .Alpha : This variable gives the mixture weights for the
%   mixture model. .Alpha(k) gives the k-th mixture weight.
%
%  - .Ea   : This variable gives the (E)xpected value of the scale
%   parameter (a) in time. This variable is denoted
%   as a in my dissertation. There is an expectation for each trajectory i
%   and each cluster k. Hence, .Ea(i,k) gives the
%   expected scale parameter for the i-th trajectory, under the k-th
%   cluster. GETCINDX(model.Ea, model.C) returns
%   the reduced n-1 matrix containing the relevant cluster-specific values.
%
%  - .Eb   : This variable gives the (E)xpected value of the translation
%   parameter (b) in time. This variable is denoted
%   as b in my dissertation. There is an expectation for each trajectory i
%   and each cluster k. Hence, .Eb(i,k) gives the
%   expected scale parameter for the i-th trajectory, under the k-th
%   cluster. GETCINDX(model.Eb, model.C) returns
%   the reduced n-1 matrix containing the relevant cluster-specific values.
%
%  - .Ee   : This variable gives the (E)xpected value of the scale
%   parameter (e) in measurement space. This variable is denoted
%   as c in my dissertation. There is an expectation for each trajectory i,
%   in each dimension d, for each cluster k. Hence, .Ee(i,k,d) gives the
%   expected scale parameter for the i-th trajectory, under the k-th
%   cluster in the d-th dimension. GETCINDX(model.Ee, model.C) returns
%   the reduced n-D matrix containing the relevant cluster-specific values.
%
%  - .Ef   : This variable gives the (E)xpected value of the translation
%   parameter (f) in measurement space. This variable is denoted
%   as d in my dissertation. The layout is the same as for .Ee.
%   GETCINDX(model.Ef, model.C) returns the reduced n-D matrix 
%   containing the relevant cluster-specific values.
%
%  - .Va   : This variable gives the (V)ariance of the scale parameter (a)
%   in time. This variable is denoted as r in my dissertation.
%   The layout is the same as for .Ea.
%   GETCINDX(model.Va, model.C) returns the reduced n-1 matrix 
%   containing the relevant cluster-specific values.
%
%  - .Vb   : This variable gives the (V)ariance of the translation
%   parameter (b) in time. This variable is denoted
%   as s in my dissertation. The layout is the same as for .Eb.
%   GETCINDX(model.Vb, model.C) returns the reduced n-1 matrix 
%   containing the relevant cluster-specific values.
%
%  - .Ve   : This variable gives the (V)ariance of the scale parameter (e)
%   in measurement space. This variable is denoted as u in my dissertation.
%   The layout is the same as for .Ee.
%   GETCINDX(model.Ve, model.C) returns the reduced n-D matrix 
%   containing the relevant cluster-specific values.
%
%  - .Vf   : This variable gives the (V)ariance of the translation
%   parameter (f) in measurement space. This variable is denoted
%   as v in my dissertation. The layout is the same as for .Ee.
%   GETCINDX(model.Vf, model.C) returns the reduced n-D matrix 
%   containing the relevant cluster-specific values.
%
%  - .Vab  : This variable gives the co-(V)ariance of (a) and (b). This
%   variable is denoted as Vab in my dissertation. The layout is the 
%   same as for .Ea.
%   GETCINDX(model.Vab, model.C) returns the reduced n-1 matrix 
%   containing the relevant cluster-specific values.
%
%  - .Vef  : This variable gives the co-(V)ariance of (e) and (f). This
%   variable is denoted as Vcd in my dissertation. The layout is the 
%   same as for .Ee.
%   GETCINDX(model.Vef, model.C) returns the reduced n-D matrix 
%   containing the relevant cluster-specific values.
%
%  - .R    : This variable gives the variance parameter for the learned
%   prior model of the scale parameter a. 
%
%  - .S    : This variable gives the variance parameter for the learned
%   prior model of the translation parameter b. 
%
%  - .T    : This variable gives the variance parameter for the learned
%   prior model of the translation parameter f. .T(k,d) gives the value for
%   the k-th cluster in the d-th dimension.
%
%  - .U    : This variable gives the variance parameter for the learned
%   prior model of the scale parameter e. .U(k,d) gives the value for
%   the k-th cluster in the d-th dimension.


PROGNAME = 'cct_model_struct';
if (exist('x')~=1)
  helpwin(PROGNAME);
else
  more on;
  help(PROGNAME)
  more off;
end

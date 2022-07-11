function cct_background(x)
% CCT_Background
%  The Curve Clustering Toolbox 1.0
%   Background on Joint Curve Clustering and Alignment
%
%   The focus of this toolbox is on the introduction of joint curve
%   clustering and alignment algorithms. Curves are often misaligned in 
%   real-world data sets. For example, the growth acceleration curves for 
%   a set of boys might share a common shape or structure due to similarities 
%   in human growth development but may be significantly misaligned due to 
%   differences in individual growth dynamics. Clustering and prediction in 
%   this situation can be difficult.
%
%   With this type of curve data, typically one of two sequential strategies 
%   is followed: (1) clustering is followed by a within-cluster alignment, 
%   or (2) the curves are first aligned globally and then clustering is carried 
%   out on the aligned curves. Both of these methods are potentially suboptimal 
%   since the clustering and alignment are often dependent on each other. The 
%   joint clustering-alignment models in this toolbox address this problem by 
%   simultaneously solving both the clustering and alignment problems using a 
%   joint EM framework.
% 
%   Standard Regression Mixture Model Specification
% 
%   The general model specification used throughout this toolbox is based
%   on the standard polynomial regression model:
%            y = XB + e,
%   where y is an observed curve, X is the associated regression matrix, B is
%   a vector of regression coefficients, and e is an error term (zero-mean
%   Gaussian). In the case of splines, X is a spline basis matrix and B is the
%   spline coefficient vector. 
%
%   If you repeat this model over K cluster components (or groups), then the
%   resulting overall density of y can be described by a linear regression 
%   mixture (LRM) model. (Some authors may choose to call this a polynomial 
%   regression mixture (PRM) model if y is a polynomial in x, but we just
%   stick with the notation LRM in this toolbox.) The mixture density for y
%   can be written as
%             p(y|Q) = \sum_k w_k N(y|Q_k),
%   where Q and Q_k denote an appropriate set of parameters for the particular
%   chosen error model (Gaussian), w_k is the k-th mixture weight, and N() is
%   the Gaussian density (note: \sum_k is just a sum indexed from 1 to K).
%   There is a straightforward EM procedure to estimate Q and determine the
%   cluster memberships for each y-curve. The EM procedure iterates
%   between estimating the cluster memberships and solving K weighted
%   least-squares problems. The algorithm for this procedure is contained in the
%   function lrm() for this toolbox. The function curve_clust() provides a
%   standard calling interface that should be used to call this and all other 
%   EM algorithms in this toolbox.
%
%   Curve-Aligned Clustering Model Specification
% 
%   In brief overview, simultaneous curve alignment and clustering is accounted
%   for by adding up to four scalar transformation variables to the above 
%   regression mixtures specification. We denote these variables by a,b,c, and d. 
%   Below we refer to x as time (and refer to transformations in time) and y as 
%   measurements (and refer to transformations in measurement space).
%   The transformed regression model can be written as
%
%            y = c[ax+b]B + d + e,
%
%   where the square brackets '[' and ']' denote the transformed regression
%   matrix (discussed below). The variable c allows for scaling in measurement
%   space, d allows for translation in measurement space, a allows for scaling
%   in time, and b allows for translation in time. This can be rewritten to more
%   clearly separate the actions of a,b and c,d:
%
%           (y-d)
%           -----  = [ax+b]B + e*,
%             c
%
%   where e* is a new error term. From this equation it is clear that c and d
%   "operate" on y (in measurement space), and a and b "operate" on x (in time).
%
%   The transformed regression matrix [ax+b] can be computed by applying the
%   transformation ax+b to the observation times before the regression matrix 
%   is built. The standard regression matrix in polynomial regression is the
%   Vandermonde matrix:
%
%           X = (1  x  x^2  x^3 ... x^p),
%
%   where x = (x_1,x_2,...,x_n)' and x^m is the vector x whose n components are
%   each raised to the power m. Thus, in this case, X is an n-by-(p+1) regression 
%   matrix. Using this notation we can easily describe [[ax+b]] as follows:
%
%        [ax+b] = (1  ax+b  (ax+b)^2  (ax+b)^3  ...  (ax+b)^p),
%
%   where (ax+b)^m is the vector ax+b whose n components are each raised to the
%   power m.
%
%   If we repeat the transformed regression model for y over K cluster components,
%   then the resulting overall density for y is a curve-aligned regression
%   mixture model. A set of linear EM algorithms (linear in the number of curve 
%   measurements) that simultaneously solves for cluster memberships and curve
%   alignments for this model was introduced in (Gaffney, 04). These joint
%   clustering-alignment EM algorithms make up the body of what the CCToolbox
%   provides. The basic set of cluster models provided in this toolbox are 
%   described in the CCT_Cluster_Models help file. See also CCT_CLUSTER_MODELS.

%  Scott J. Gaffney   18 May 2005
%  University of California, Irvine

PROGNAME = 'cct_background';
if (exist('x')~=1)
  helpwin(PROGNAME);
else
  more on;
  help(PROGNAME)
  more off;
end


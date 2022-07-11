function cct_curve_format(x)
%CCT_Curve_Format 
% Curve Clustering Toolbox 1.0
%  Supported Data Structure Formats
%
%  Curve representation
%  ====================
%  Curves are represented throughout the CCToolbox in several different
%  ways. Some formats are best for computation whereas others are best for
%  ease-of-use. The three main formats (Cell, Feature Vector, and Sequence)
%  used in this toolbox are discussed below; however, the standard user need
%  only be concerned with the Cell format and how to pass curve data to
%  functions using matlab structures.
%
%  The 'Trajs' structure
%  ---------------------
%    All user functions in the CCToolbox accept curve data through the standard
%  argument named 'Trajs'. In common use, one should pass-in a structure for
%  this argument with the following two fields:
%    Trajs
%     .Y   : observed curve values
%     .X   : observation time-points at which Trajs.Y was observed
%  The field 'X' is optional in which case a default value for it is assumed as
%  described below. It is also possible to just pass-in Trajs.Y directly.
%  In other words, you can dispense with the structure altogether if you are
%  only interested in the curve values and not the time-points themselves. The
%  user functions in the CCToolbox will correctly detect whether or not a 
%  structure has been used.
%    The 'X' and 'Y' fields can be given in either the Cell or Feature Vector
%  formats described below; the Cell format is preferred. Before we detail the
%  curve formats, we list some general curve notation.
%    Assume that we have curves in the form y=f(x) so that the observed curve
%  values y are functions of an independent variable x (usually associated
%  with time). The curve values y can be multivariate but x is assumed to
%  be univariate.
%    n  : the number of curves
%    ni : length of the i-th curve
%    D  : the dimension of y (e.g. D is 1 if y is univariate)
%    N  : the total number of "points" in all curves, this is equal to the
%         sum of the lengths of all curves (independent of the dimension of y)
%
%  The Cell Format
%  ---------------
%  In this format curves are represented using cell arrays. Let Y be an 
%  n-by-1 cell array. Then Y{i} is an ni-by-D matrix holding the observed curve
%  values for the i-th curve.  If the observation times (i.e., x) are uniform 
%  and regular across all curves, then there is no need to specify
%  or keep track of X directly. In this case, each curve is assumed to have
%  time points [0,1,...,(ni-1)]. Otherwise, X is an n-by-1 cell array and X{i} is
%  a vector of length ni holding the time points for the i-th curve. 
%    As a special case, if X is a single vector, then it is assumed that 
%  all curves share the specified time points up to their length ni. In 
%  this case, the length of the vector X must be greater than or equal to 
%  the longest curve in Y. 
%
%  The Feature Vector Format
%  -------------------------
%  Sets of curves are represented by a single matrix. The n curves
%  are all assumed to have uniform length m and dimension D. Thus the curve
%  values can be contained in the matrix Y of size n-by-m-by-D. The observation
%  times are not usually important when curves are represented in this format.
%  However, the toolbox will assume the time points [0,1,...,(m-1)] if needed.
%  Similarly, X can be a single vector as described above in Section "The Cell
%  Format".
%
%  The Sequence Format
%  -------------------
%  This format represents curves using three matrices: Y, X, and Seq.
%  This format is used internally for cluster algorithms in the CCToolbox 
%  and is not usually embedded in a 'Trajs' structure for user functions.
%  The three matrices are laid out as follows:
%    (1) Y is an N-by-D matrix holding the N observed curve values (each of
%        dimension D), one curve after another.
%    (2) X is an N-by-1 matrix holding the N time-points at which the curve 
%        values in Y were observed. In advanced usage, X can also be 
%        an N-by-(p+1) regression matrix where p gives the order of the 
%        regression (i.e., the matrix will be used directly without 
%        modification). However, this can only be used when no time
%        transformations are allowed (i.e., lrm, srm, lrm_d, srm_d, lrm_cd,
%        srm_cd.)
%    (3) Seq is a vector of length (n+1) holding the curve offsets into X 
%        and Y (e.g., Seq(3) gives the starting index of the third curve in 
%        both X and Y). Seq(n+1) is ALWAYS equal to Seq(n) + nj, where
%        nj equals the length of the last curve. In other words, Seq(n+1)
%        points to the position just beyond the last curve entry. This
%        means that Seq(n+1)-1 is always equal to N (the total number of
%        points, as defined above).
%
%
%  Using the different curve formats
%  ---------------------------------
%  As far as the common user is concerned, there is no need to use anything
%  other than the standard 'Trajs' structure and the Cell format.
%  However, some functions only operate on or return specific formats. 
%  One can usually tell what formats are acceptable by noting the name of 
%  the input/output arguments for a function.
%    Functions receive or return curve data using two different parameter sets:
%  (1) either the single argument 'Trajs', or (2) the three arguments Y,
%  X, and Seq. In the latter case, only the sequence format is acceptable. In
%  the former case (when 'Trajs' is specified), the function operates as
%  described above.


PROGNAME = 'cct_curve_format';
if (exist('x')~=1)
  helpwin(PROGNAME);
else
  more on;
  help(PROGNAME)
  more off;
end

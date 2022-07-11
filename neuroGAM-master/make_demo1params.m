% funtion definitions
function prs = make_demo1params(sr)
%% define parameters
prs=[];
prs.varname =  {{'Position-x (m)' , 'Position-y (m)'} ,{'View-x (m)' , 'View-y (m)' }, 'Position-z (m)' ,'View-z (m)', 'Angular speed (deg/s)' , 'Body Speed (cm/s)'};
prs.vartype = {'2D' '2D' '1D' '1D' '1D' '1D'};
prs.basistype={'raisedcosine','raisedcosine','raisedcosine','raisedcosine','raisedcosine','raisedcosine'};
prs.nbins = {[20 , 20],[20 , 20] ,6, 6, 10 , 10};
prs.binrange = [];%{[0 0.6 ;0 1.2],[0 0.6 ;0 1.2] ,[0.2; 1.4],[0.2 ;1.4],[0;2000],[0 ;300] };
prs.nfolds = 5;
prs.dt = 1/sr;
prs.filtwidth = 12;
prs.linkfunc = 'log';
prs.lambda = {10 ,10, 10, 10, 50 ,50};
prs.alpha = 0.05;
prs.varchoose=[0 0 0 0 0 0];
prs.method='FastForward';
end
%Problems to solve with this GAM toolbox:
%Describe types of basistype input
%Fix conv3 function requirement
%for 2d data bins size product should be a square (equal sized bins)
%I had to comment out from line 153 in BuildGAM.m
%output of goodness of fit is  [varExplain_test varExplainPseudo_test correlation_test log_llh_test mse_test sum(n) numel(test_ind)];
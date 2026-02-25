function R=faster_corr_mtrx(data)
%function R=faster_corr_mtrx(data)
%
% A much faster way to compute a correlation matrix than MATLAB's built in
% corr.m function. 
%
% Input:
%  data - Column vector or matrix of data (columns are observations, rows
%         are variables)
%
% Output:
%  R - Matrix of linear correlation coefficients (the xth, yth element is
%      the linear correlation coeffecient between variables x and y).
%
% Author:
% David Groppe
%

C=cov(data);
n_dim=length(C);
R=zeros(n_dim,n_dim);
sd=sqrt(diag(C));
for a=1:n_dim,
    for b=1:n_dim,
        R(a,b)=C(a,b)/(sd(a)*sd(b));
    end
end
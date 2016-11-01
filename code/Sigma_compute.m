function [Sigma, Sigma_svd] = Sigma_compute(X)
% Compute the sample nonparanormal covariance matrice for the two groups of
% subjects and the K views of data. 
%
%   [Sigma, Sigma_svd] = Sigma_compute(X) compue the  the sample
%   nonparanormal covariance matrice using the gene expression profiels X.
%
%   INPUT:
%       X: K-by-2 cell, where each element, X{kc}, is a n_c-by-p matrix that
%       represents the gene expression profie for the k-th view of data and
%       the c-th group of subject. 
%
%   Output:
%       Sigma: K-by-2 cell, where each element, Sigma{kc}, is a p-by-p matrix that
%       represents the sample nonparanormal covariance matrix for the k-th
%       view of data and the c-th group of subject.
%       Sigma_svd: a cell with length K which provides the svd
%       decomposition of  A = (Sigma{k,1}'*Sigma{k,2} +
%       Sigma{k,2}'*Sigma{k,1})/2. That is,  [Sigma_svd{k}.U, 
%       Sigma_svd{k}.D] = svd(A), and  Sigma_svd{k}.D =
%       diag(Sigma_svd{k}.D)
%
% See also pDNA, pDNA_admm_iters, solve_G, soft_thresh
%
%
% Reference
%   X. F. Zhang, L. Ou-Yang, and H yan (2016) Incorporating prior
%   information into differential network analysis using graphical models
%
% COPYRIGHT  Central China Normal University
% Xiao-Fei Zhang <zhangxf@mail.ccnu.edu.cn>   

%%
% compute the number of views and the number of genes
K = size(X,1);
p = size(X{1,1},2);

Sigma = cell(K,2);
Sigma_svd = cell(K,1);

for k = 1:K
   
    
    % compute the kendall tau correlation using Equation (1) in the main
    % text
    Sigma{k,1} =  kendalltau(X{k,1});
    Sigma{k,2} =  kendalltau(X{k,2});
    Sigma{k,1}(isnan(Sigma{k,1})) = 0;
    Sigma{k,2}(isnan(Sigma{k,2})) = 0;
    
    % the sample nonparanormal covariance matrix using Equation (2) in the
    % main text
    Sigma{k,1} = sin(0.5*pi*Sigma{k,1});
    Sigma{k,2} = sin(0.5*pi*Sigma{k,2});
    
    % set the diagonal element to be 1
    Sigma{k,1} =  Sigma{k,1} - diag(diag(Sigma{k,1})) + eye(p);
    Sigma{k,2} =  Sigma{k,2} - diag(diag(Sigma{k,2})) + eye(p);
    
    % compute the svd decomposition of A
    A = (Sigma{k,1}'*Sigma{k,2} + Sigma{k,2}'*Sigma{k,1})/2;
    [Sigma_svd{k}.U, Sigma_svd{k}.D] = svd(A);
    Sigma_svd{k}.D = diag(Sigma_svd{k}.D);
end

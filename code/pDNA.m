function [Delta, V] = pDNA(Sigma, lambda, varargin)

% Complete algorithm of the pDNA model which is provided in Algorithm 1 in
% the Supplementary Data.
%
%   [Delta, V] = pDNA(Sigma, lambda) solves the pDNA model using the sample
%   nonparanormal covariance matrices Sigma and the penalty paramter
%   lambda. The result Delta is the estimated precision matrix difference.
%   V is the decomposition of Delta, that is, Delta{k} = V{k} + V{k}'.
%
%  [Delta, V] = pDNA(Sigma, lambda, 'PARAM1', val1, 'PARAM2', val2...)
%  allows you to specify optional parameters to control the model solution.
%  Available parameter name/value pairs are:
%
%       'Sigma_svd': a cell with length K which provides the svd
%       decomposition of  A = (Sigma{k,1}'*Sigma{k,2} +
%       Sigma{k,2}'*Sigma{k,1})/2. That is,  [Sigma_svd{k}.U,
%       Sigma_svd{k}.D] = svd(A), and  Sigma_svd{k}.D =
%       diag(Sigma_svd{k}.D), defualt is empty
%       'F': co-pathway indication matrix,  default is a one
%       matrix
%       'out_MAX_ITER': maximum iteration of outer iterations, default is
%       50
%       'rho': the intitial value of the penalty paramter for the ADMM
%       algorithm, default is 1e-1
%       'rhoIncr': the increase size of rho, default is 1.1
%       'rhoMax': the maximize value of rho, default is 1e10
%       'MAX_ITER': maximum iteration of the ADMM algorithm, default is 5e2
%       'RELTOL': relative tolerence, default is 1e-4
%
%
%
%   INPUT:
%       Sigma: K-by-2 cell, where each element, Sigma{kc}, is a p-by-p matrix that
%       represents the sample nonparanormal covariance matrix for the k-th
%       view of data and the c-th group of subject.
%       lambda: nonnegative tuning parameter that controls the sparsity of
%       the resulting differential networks.
%
%   Output:
%       Delta: K-by-1 cell, where each element, Delta{k}, is a p-by-p
%       matrix that represents the estimated precision matrix difference
%       for the k-th view of data.
%       V: K-by-1 cell, where each element, V{k}, is a p-by-p
%       matrix that represents the decomposition of Delta{k}, that is
%       Delta{k} = V{k} + V{k}'.
%
%
% Examples
%
% See also pDNA_admm_iters, solve_G, soft_thresh, Sigma_compute
%
%
% Reference
%   X. F. Zhang, L. Ou-Yang, and H yan (2016) Incorporating prior
%   information into differential network analysis using graphical models
%
% COPYRIGHT  Central China Normal University
% Xiao-Fei Zhang <zhangxf@mail.ccnu.edu.cn>

%%
% parse inputs

argin = inputParser;
argin.addRequired('Sigma', @(x)  iscell(x));
argin.addRequired('lambda', @(x) isnumeric(x) && x>=0);

argin.addParamValue('Sigma_svd', [], @(x) iscell(x));
argin.addParamValue('F', [], @(x) isnumeric(x));

argin.addParamValue('out_MAX_ITER', 50, @(x) isnumeric(x) && x>0);
argin.addParamValue('rho', 1e-1, @(x) isnumeric(x) && x>0);
argin.addParamValue('rhoIncr', 1.1, @(x) isnumeric(x) && x>0);
argin.addParamValue('rhoMax', 1e10, @(x) isnumeric(x) && x>0);
argin.addParamValue('MAX_ITER', 5e2, @(x) isnumeric(x) && x>0);
argin.addParamValue('RELTOL', 1e-4, @(x) isnumeric(x) && x>0);
argin.parse(Sigma, lambda, varargin{:});

% Copy from params object
Sigma_svd = argin.Results.Sigma_svd;
F = argin.Results.F;
out_MAX_ITER = argin.Results.out_MAX_ITER;
rho = argin.Results.rho;
rhoIncr = argin.Results.rhoIncr;
rhoMax = argin.Results.rhoMax;
MAX_ITER = argin.Results.MAX_ITER;
RELTOL = argin.Results.RELTOL;

%%

% compute the number of views K and the number of genes p
K = size(Sigma,1);
p = size(Sigma{1,1},1);

% if the input of svd decomposistion is empty, compute the svd
% decomposistion
if isempty(Sigma_svd)
    Sigma_svd  = cell(K,1);
    for k = 1:K
        A = (Sigma{k,1}'*Sigma{k,2} + Sigma{k,2}'*Sigma{k,1})/2;
        [Sigma_svd{k}.U, Sigma_svd{k}.D] = svd(A);
        Sigma_svd{k}.D = diag(Sigma_svd{k}.D);
    end
end

% if there is no co-pathway indication matrix available, set F = ones(p,p).
% That is, we do not consider the pathway-based constains if we cannot
% obatin the pathway information.
if isempty(F)
    F = ones(p,p);
end



%%
% initial estimate by setting the weights equals to 1.
Delta = cell(K,1);
V = cell(K,1);
psi = ones(p,1);
omega  = ones(p,p);
for k = 1:K
    [Delta{k}, V{k}]  = pDNA_admm_iters(Sigma{k,1}, Sigma{k,2}, F, lambda*4, psi, omega, 'Sigma_svd', Sigma_svd{k}, 'rho', rho, 'rhoIncr', rhoIncr, ...
        'rhoMax', rhoMax, 'MAX_ITER', MAX_ITER, 'RELTOL', RELTOL);
end


%
%%
for i = 1:out_MAX_ITER
    i
    Delta_old = Delta;
    
    % compute weights psi and omega
    temp = zeros(p,p);
    for k = 1:K
        temp =  temp +  abs(V{k});
    end
    temp = sqrt(temp);
    temp1 = sqrt(sum(temp,1)');
    
    omega = 1 ./ (temp + 1e-4);
    psi = 1 ./ (temp1 + 1e-4);
    
    % solve proble (3) using ADMM algorithm
    for k = 1:K
        [Delta{k}, V{k}]  = pDNA_admm_iters(Sigma{k,1}, Sigma{k,2}, F, lambda, psi, omega, 'Sigma_svd', Sigma_svd{k}, 'rho', rho, 'rhoIncr', rhoIncr, ...
            'rhoMax', rhoMax, 'MAX_ITER', MAX_ITER, 'RELTOL', RELTOL);
    end
    
    % check the convergence condition
    R_err = 0;
    R_norm = 0;
    for k = 1:K
        R_err = R_err + norm(Delta{k} - Delta_old{k});
        R_norm = R_norm + norm(Delta_old{k});
    end  
    if R_err < R_norm * RELTOL
        break
    end
end

% output Delta 
for k = 1:K
    Delta{k} = V{k} + V{k}';
end



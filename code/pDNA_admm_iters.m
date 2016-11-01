function [Delta, V, history] = pDNA_admm_iters(Sigma_1, Sigma_2, F, lambda, psi, omega, varargin)

% ADMM algorithm for the optimization problem (3) which is provided in
% Algorithm 2 in the Supplementary Data.
%
%  [Delta, V, history] = DNA_admm_iters(Sigma_1, Sigma_2, F, lambda, psi,
%   omega) solves problem (3).
%
%  [Delta, V, history] = DNA_admm_iters(Sigma_1, Sigma_2, F, lambda, psi,
%   omega, 'PARAM1', val1, 'PARAM2', val2...) allows you to specify
%   optional parameters to control the model solution. Available parameter
%   name/value pairs are: 
%
%       'Sigma_svd': a cell with length K which provides the svd
%       decomposition of  A = (Sigma{k,1}'*Sigma{k,2} +
%       Sigma{k,2}'*Sigma{k,1})/2. That is,  [Sigma_svd{k}.U,
%       Sigma_svd{k}.D] = svd(A), and  Sigma_svd{k}.D =
%       diag(Sigma_svd{k}.D), defualt is empty
%       'Delta_init':  the intitial value of Delta,  default is empty
%       'V_init':  the intitial value of V,  default is empty
%       'rho': the intitial value of the penalty paramter for the ADMM
%       algorithm, default is 1e-1
%       'rhoIncr': the increase size of rho, default is 1.1
%       'rhoMax': the maximize value of rho, default is 1e10
%       'MAX_ITER': maximum iteration of the ADMM algorithm, default is 5e2
%       'RELTOL': relative tolerence, default is 1e-4
%
%
%   INPUT:
%       Sigma_1: a p-by-p matrix that represents the sample nonparanormal
%       covariance matrix for the first group of subjects
%       Sigma_2: a p-by-p matrix that represents the sample nonparanormal
%       covariance matrix for the second group of subjects
%       F: a p-by-p matrix that represents the co-pathway indication
%       matrix
%       lambda: nonnegative tuning parameter that controls the sparsity of
%       the resulting differential networks.
%       psi: a p-by-1 vector that represents the weights.
%       omega: a p-by-p matrix that represents the weights.
%
%
%   Output:
%       Delta: a p-by-p matrix that represents the estimated precision
%       matrix difference
%       V: a p-by-p matrix that represents the decomposition of Delta, that is
%       Delta = V + V'.
%
%
% Examples
%
% See also pDNA, solve_G, soft_thresh, Sigma_compute
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
argin.addRequired('Sigma_1', @(x)  isnumeric(x));
argin.addRequired('Sigma_2', @(x)  isnumeric(x));
argin.addRequired('F', @(x)  isnumeric(x));
argin.addRequired('lambda', @(x) isnumeric(x)&& x>=0 );
argin.addRequired('psi', @(x) isnumeric(x) );
argin.addRequired('omega', @(x) isnumeric(x) );


argin.addParamValue('Sigma_svd', [], @(x) isstruct(x));
argin.addParamValue('Delta_init', [], @(x) isnumeric(x));
argin.addParamValue('V_init', [], @(x) isnumeric(x));
argin.addParamValue('rho', 1e-1, @(x) isnumeric(x) && x>0);
argin.addParamValue('rhoIncr', 1.1, @(x) isnumeric(x) && x>0);
argin.addParamValue('rhoMax', 1e10, @(x) isnumeric(x) && x>0);
argin.addParamValue('MAX_ITER', 5e2, @(x) isnumeric(x) && x>0);
argin.addParamValue('RELTOL', 1e-4, @(x) isnumeric(x) && x>0);
argin.parse(Sigma_1, Sigma_2, F, lambda, psi, omega, varargin{:});

% Copy from params object
Sigma_svd = argin.Results.Sigma_svd;
Delta_init = argin.Results.Delta_init;
V_init = argin.Results.V_init;
rho = argin.Results.rho;
rhoIncr = argin.Results.rhoIncr;
rhoMax = argin.Results.rhoMax;
MAX_ITER = argin.Results.MAX_ITER;
RELTOL = argin.Results.RELTOL;

%% initialize

% compute the number of genes
p = size(Sigma_1,1);

% if the input of svd decomposistion is empty, compute the svd
% decomposistion
if isempty(Sigma_svd)
    A = (Sigma_1'*Sigma_2 + Sigma_2'*Sigma_1)/2;
    [Sigma_svd.U, Sigma_svd.D] = svd(A);
    Sigma_svd.D = diag(Sigma_svd.D);
end


% initialize Delta
if isempty(Delta_init)
    Delta = eye(p);
else
    Delta = (Delta_init + Delta_init')/2;
end

% initialize V
if isempty(V_init)
    V = eye(p);
else
    V = V_init;
end

% initialize W
W = V';

% initialize dual variables to be zero matrices
P = zeros(p,p);
Q = zeros(p,p);

% rewrite psi as a matrix form
psi = repmat(psi',p, 1);

%% ADMM iteration

for iter = 1:MAX_ITER
    iter;
    Delta_old = Delta;
    
    % updata Delta accoring to Equation (8) in the supplementary data
    Delta = solve_G(Sigma_svd.U, Sigma_svd.D, Sigma_1 - Sigma_2 + rho*(V+W) - P, rho);
    
    % updata V accoring to Equation (11) in the supplementary data
    H = 0.5*(Delta - W + W') + 0.5*(P - Q)/rho;
    V = (soft_thresh(H, lambda.*psi.*omega/(8*rho))).*F;
    
    % updata W accoring to Equation (12) in the supplementary data
    W = 0.5*(Delta - V + V') + 0.5*(P + Q')/rho;
    
    % updata dual variables using a dual-ascend updata rule
    P = P + rho*(Delta - (V + W));
    Q = Q + rho*(V - W');
    
    % check the convergence condition according to Equation (13-15) in the supplementary data
    history.rho(iter) = rho;
    history.Delta_pri_diff(iter)  = sum(sum( abs(Delta - Delta_old))) ./ (sum(sum(abs(Delta_old)))+eps);
    history.Delta_dual_diff(iter) = sum(sum( abs(Delta - (V + W)))) ./ (sum(sum( abs(Delta)))+eps);
    history.V_dual_diff(iter) = sum(sum( abs(V - W'))) ./(sum(sum( abs(V)))+eps);
    %     history.Delta_pri_diff(iter)  = norm(Delta - Delta_old) ./ (norm(Delta_old)+eps);
    %     history.Delta_dual_diff(iter) = norm(Delta - (V + W)) ./ (norm(Delta_old)+eps);
    %     history.V_dual_diff(iter) = norm(V - W') ./ (norm(V)+eps);
    if max([ history.Delta_pri_diff(iter), history.Delta_dual_diff(iter), history.V_dual_diff(iter)]) < RELTOL
        break;
    end
    
    % update rho
    rho = min(rho*rhoIncr, rhoMax);
end


end

function X = solve_G(UA, DA, B, rho)
% Operator G(A, B) is defined in Equation (6) in the Supplementary
% Data. Here A = UA*diag(DA)*UA', and this function is used to compute G(A
% + rho*I, B). It is well known that A + rho*I =  UA * (diag(DA)+  rho*I) *
% UA' 
% 

% See also pDNA, pDNA_admm_iters, soft_thresh, Sigma_compute
%
%
% Reference
%   X. F. Zhang, L. Ou-Yang, and H yan (2016) Incorporating prior
%   information into differential network analysis using graphical models
%
% COPYRIGHT  Central China Normal University
% Xiao-Fei Zhang <zhangxf@mail.ccnu.edu.cn>


p = length(DA);
DA_rho = DA + rho;

C = 2./(repmat(DA_rho,1,p)+(repmat(DA_rho,1,p))'); 
X = UA*((UA'*B*UA).*C)*UA';


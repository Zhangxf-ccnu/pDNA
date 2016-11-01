%% -------------------- SOFT-THRESHOLDING OPERATOR ------------------- %%
%
%
% ---------------- minimize 1/2*||X - Y||_2^2 + lambda*||X||_1 -------- %
%
% See also pDNA, pDNA_admm_iters, solve_G, Sigma_compute
%
%
% Reference
%   X. F. Zhang, L. Ou-Yang, and H yan (2016) Incorporating prior
%   information into differential network analysis using graphical models
%
% COPYRIGHT  Central China Normal University
% Xiao-Fei Zhang <zhangxf@mail.ccnu.edu.cn>   


function [X] = soft_thresh(Y,lambda)

if nargin == 1
    disp('Not enough inputs');
    disp('Enter both Y and lambda');
    return
end
        
X = sign(Y).*(max(abs(Y) - lambda,0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			TEST FILE FOR pDNA using glioblastoma gene expression data                     %
%  Refer to the paper: X. F. Zhang, L. Ou-Yang, and H yan (2016)
%  Incorporating prior information into differential network analysis using nonparanormal graphical models  %
% for further details on the pDNA model.                      %	
%                                                                                       %
% CONTACT Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) for any questions or comments on the code.			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


% load data
load GBM_Exp
lambda  = 0.55;


% estimate the differential networks using pDNA (with lambda = 0.55)
[Delta_hat, V_hat] = pDNA(GBM_Exp.Sigma, lambda, 'Sigma_svd', GBM_Exp.Sigma_svd, 'F', GBM_Exp.Co_pathway);

%%
p = size(Delta_hat{1},1);
K = length(Delta_hat);

% compue the weighted differential network
W = zeros(p,p);
for k = 1:K
    W = W + double(Delta_hat{k}~=0);
end

% compute the degrees of nodes
Degree = sum(W,2);

[~,ID] = sort(Degree,'descend');

% show the weighted differential network accoring to the degrees of nodes
W = W(ID, ID);
imagesc(W)
Gene_sorted = GBM_Exp.Gene(ID);

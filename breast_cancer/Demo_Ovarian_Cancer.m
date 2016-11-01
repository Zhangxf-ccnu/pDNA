%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			TEST FILE FOR pDNA using ovarian gene expression data                     %
%  Refer to the paper: X. F. Zhang, L. Ou-Yang, and H yan (2016)
%  Incorporating prior information into differential network analysis using nonparanormal graphical models  %
% for further details on the pDNA model.                      %	
%                                                                                       %
% CONTACT Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) for any questions or comments on the code.			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


% load data
load OV_Exp
lambda  = 0.35;


% estimate the differential networks using pDNA (with lambda = 0.35)
[Delta_hat, V_hat] = pDNA(OV_Exp.Sigma, lambda, 'Sigma_svd', OV_Exp.Sigma_svd, 'F', OV_Exp.Co_pathway);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 			TEST FILE FOR pDNA using simulated data                     %
%  Refer to the paper: X. F. Zhang, L. Ou-Yang, and H yan (2016)
%  Incorporating prior information into differential network analysis using graphical models  %
% for further details on the pDNA model.                      %	
%                                                                                       %
% CONTACT Xiao-Fei Zhang (zhangxf@mail.ccnu.edu.cn) for any questions or comments on the code.			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% experiment settings
clear 
close all 

K = 6;
p = 100;
n_1 = 200;
n_2 = 200;
m_pert = 4;
rho_1 = 0.5;
rho_2 = 0.9;
umin_sparse = 0.5;
umax_sparse = 1;

%% generate simulation data
[X, Delta_true, F] = generate_data(K, p, n_1, n_2, m_pert, rho_1, rho_2,  umin_sparse, umax_sparse);

%% compute the  sample nonparanormal covariance matrice
[Sigma, Sigma_svd] = Sigma_compute(X);

%% estimate the differential netwoks using pDNA
[Delta, V] = pDNA(Sigma, 0.8,'F',F);

%% Compare estimated networks to true networks

subplot(2,3,1);
imagesc(Delta_true{1} - diag(diag(Delta_true{1})));
colorbar
title('True \Delta^1');


subplot(2,3,2);
imagesc(Delta{1}-diag(diag(Delta{1})))
colorbar
title('Estimated \Delta^1');


subplot(2,3,3);
imagesc(V{1}-diag(diag(V{1})))
colorbar
title('Estimated V^{1}');


subplot(2,3,4);
imagesc(Delta_true{2} - diag(diag(Delta_true{2})));
colorbar
title('True \Delta^2');

subplot(2,3,5);
imagesc(Delta{2}-diag(diag(Delta{2})))
colorbar
title('Estimated \Delta^2');

subplot(2,3,6);
imagesc(V{2}-diag(diag(V{2})))
colorbar
title('Estimated V^{2}');

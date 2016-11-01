function [X, Delta, F] = generate_data(K, p, n_1, n_2, m_pert, rho_1, rho_2,  umin_sparse, umax_sparse)
% Generate the simulation data using the procedure provised in Algorithm 3
% in the Supplementary Data. For details of the input and output, please
% refer to  Algorithm 3.
%
% Examples
%
% See also pDNA, pDNA_admm_iters, solve_G, soft_thresh, Sigma_compute
%
%
% Reference
%   X. F. Zhang, L. Ou-Yang, and H yan (2016) Incorporating prior information into differential network analysis using nonparanormal graphical models.
%
% COPYRIGHT  Central China Normal University
% Xiao-Fei Zhang <zhangxf@mail.ccnu.edu.cn>
%% step 1
Net1 = SFNG(floor(0.6*p), 2, 1);
Net2 = SFNG(floor(0.6*p), 2, 1);
Net = zeros(p,p);
Net(1:floor(0.6*p),1:floor(0.6*p)) = Net1;
Net((floor(0.4*p)+1):p,(floor(0.4*p)+1):p) = Net((floor(0.4*p)+1):p,(floor(0.4*p)+1):p) + Net2;
Net(Net > 0) = 1;

%% step 2
Pathway = zeros(p,2);
Pathway(1:floor(0.6*p),1) = 1;
Pathway((floor(0.4*p)+1):p,2) = 1;
F = Pathway*Pathway';
F(F~=0) = 1;

%% step 3
ID_pert = sort(randsample(p, m_pert)) ;

ID_comm = cell(length(ID_pert),1);
for i = 1:length(ID_pert);
    ID_comm{i} = randsample(p, floor(p*rho_1*(1-rho_2)));
end

%% step 4


W = cell(K,2);
Theta_true = cell(K,2);
Delta = cell(K,1);
X = cell(K,2);


for k = 1:K
    
    % step 4(a)
    W{k,1} = (ones(p,p) * umin_sparse + rand(p,p) * (umax_sparse - umin_sparse)) .* ( 2*unidrnd(2,p,p)-3);
    W{k,1} = triu(W{k,1}, 1) + (triu(W{k,1}, 1))';
    W{k,1} =  W{k,1}.*Net;
    
    % step 4(b)
    W{k,2} =  W{k,1};
    
    % step 4(c)
    for i = 1:length(ID_pert);
        
        % step 4(c)i
        ID_specific  = randsample(setdiff(1:p, ID_comm{i}), ceil(p*rho_1*rho_2));
        % step 4(c)ii
        ID = union(ID_comm{i}, ID_specific);
        % step 4(c)iii
        temp =  ((ones(length(ID),1) * umin_sparse) + (rand(length(ID),1) * (umax_sparse - umin_sparse))) .* (2*unidrnd(2,length(ID),1)-3);
        % step 4(c)iv
        flipper = unidrnd(2);
        % step 4(c)v
        if flipper == 1
            W{k,1}(ID_pert(i),ID) = temp';
            W{k,1}(ID,ID_pert(i)) = temp;
        else
            W{k,2}(ID_pert(i),ID) = temp';
            W{k,2}(ID,ID_pert(i)) = temp;
        end
        
    end
    
    
    % step 4(d)
    W{k,1} = W{k,1}.*F;
    W{k,2} = W{k,2}.*F;
    
    % step 4(e)
    min_eigval_1 = min(eig(W{k,1}));
    min_eigval_2 = min(eig(W{k,2}));
    min_eigval = min(min_eigval_1, min_eigval_2);
    Theta_true{k,1}  = W{k,1} + (eye(p)*(abs(min_eigval) + 0.1));
    Theta_true{k,2}  = W{k,2} + (eye(p)*(abs(min_eigval) + 0.1));
    
    % step 4(f)
    Delta{k} =  Theta_true{k,2} -  Theta_true{k,1};
    
    % step 4(g)
    X{k,1} = mvnrnd(zeros(p,1),inv(Theta_true{k,1}), n_1);
    X{k,2} = mvnrnd(zeros(p,1),inv(Theta_true{k,2}), n_2);
    
end



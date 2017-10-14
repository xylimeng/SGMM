% function to return initials 
% input: y - obs (n by d), p - template (n by S) 
% if 'p' is provided : use matrix operations; otherwise, use GMM 
function par_ini = obs2ini(y, p)

if nargin == 2
    n = size(y, 1);
    S = size(p, 2); % S = number of regions; here we have three regions
    d = size(y, 2); % d = dimension of data; here we have BL and FU
    
    U = zeros(d,S); % store mu's - each col is a region
    V = zeros(d,d,S); % store covariance matrix
    for s = 1:S
        ps = p(:,s);
        % Mean by region
        U(:,s) = y' * ps / sum(ps,1);
        % Variance by region
        yc = y - repmat(U(:,s)',n,1);
        V(1,1,s) = (yc(:,1).^2)' * ps / sum(ps,1);
        V(2,2,s) = (yc(:,2).^2)' * ps / sum(ps,1);
        V(1,2,s) = (yc(:,1).*yc(:,2))' * ps / sum(ps,1);
        V(2,1,s) = V(1,2,s);
        V(:,:,s) = V(:,:,s);
    end
    
    % another way to calculate V - numerically the same
    V2 = zeros(d,d,S);
    for s = 1:S
        ps = p(:,s);
        % Mean by region
        U(:,s) = y' * ps / sum(ps,1);
        % Variance by region
        yc = y - repmat(U(:,s)',n,1);
        V2(:,:,s) = yc' * (yc .* repmat(ps, [1, 2])) /sum(ps, 1);
    end
    
    
    par_ini.mu = U;
    par_ini.sigma = V; % covariance matrix
    par_ini.p = mean(p, 1);
else
    options = statset('MaxIter', 1000);
    obj_GM = gmdistribution.fit(obs, K, 'Replicates', 2, 'Options', options);
    % extract estimates
    par_GM.mu = obj_GM.mu';
    par_GM.sigma = obj_GM.Sigma;
    par_GM.p = obj_GM.PComponents;
    par_ini = par_GM;
end

par_ini.gamma = [1/3, 1/3, 1/3]; % each col is a region: 1 by K
end

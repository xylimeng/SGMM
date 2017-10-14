% function to transfer b_{ik} to p_{ik}
% gamma is a vector with length K
% P is the inputed 2-D matrix (n by K)
function [P_transform] = P2Pt(P, gamma0)

K = size(P, 2); % number of clusters 
if nargin < 2
    gamma0 = ones([1,K]); % no transformation at all
end

P_transform = zeros(size(P)); 

K = size(P, 2);

for k = 1:K
    P_transform(:, k) = P(:,k) .* gamma0(k);
end

temp_sum = sum(P_transform,2);

for k = 1:K
    P_transform(:,k) = P_transform(:,k)./temp_sum;
end
end



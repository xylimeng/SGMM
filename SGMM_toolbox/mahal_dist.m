function mahal = mahal_dist(obs, MU, SIGMA)
% function to calculate the mahal distance - root of (obs - MU)' *
% inv(SIGMA) * (obs - MU)
% obs - n by K; MU - 1 by K; SIGMA - K by K
% mahal - n by 1
m = reshape(MU, [1, length(MU)]); 
M = m(ones(size(obs, 1),1),:);
root_C = sqrtm(inv(SIGMA)); 
a = (obs - M) * root_C; 
mahal = sqrt(sum(a .* a, 2)); 
end 

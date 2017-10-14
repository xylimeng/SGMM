% function to estimate parameters based on the segmentation results; 
% type:   hard - use labels to estimate parameters; 
%         soft - use probability maps to estimate parameters; 
%                note here we typically should use the posterior probability maps invovling the data 

% Version: 0.0 created at 23-Nov-2013

function par = EstimatePars(input, obs, type)
% EstimatePars estimates (mu, sigma) from 'input' (either estimated lables
% or the probability maps)
% the observation 'obs'.

if nargin < 3
    type = 'hard';
end

% fprintf([type, ' segmention to estimate (mu, sigma) \n'])

if strcmp(type, 'hard')
    EstLabels = input;
    estimate_labels = unique(EstLabels); % typically 1:K
    K = length(estimate_labels); % number of labels
    par.mu = zeros([1,K]);
    par.sigma = zeros([1,K]);
    par.p = zeros([1,K]);
    for i = 1:K
        par.mu(i) = mean(obs(EstLabels == estimate_labels(i)));
        par.sigma(i) = std(obs(EstLabels == estimate_labels(i)));
        par.p(i) = mean(EstLabels == estimate_labels(i));
    end
end

if strcmp(type, 'soft')
    ProbabilityMap = input;
    K = size(ProbabilityMap, 2);
    
    par.mu = zeros([1,K]);
    par.sigma = zeros([1,K]);
    par.p = zeros([1,K]);    
    for j = 1:K
        par.p(j) = mean(ProbabilityMap(:,j));
        par.mu(j) = ProbabilityMap(:,j)' * obs / sum(ProbabilityMap(:,j));
        par.sigma(j) = sqrt((ProbabilityMap(:,j)' * ((obs - par.mu(j)).^2)) / sum(ProbabilityMap(:,j)));
    end
end

end

    

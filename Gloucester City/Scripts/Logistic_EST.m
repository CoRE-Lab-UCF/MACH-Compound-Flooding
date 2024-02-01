%This estimate the estimate and CI values for the given values in XI using
%sample data 

% 
function[Xi,MLE,CI]= Logistic_EST(data,Xi,numBootstrap)
    % Check input
    if nargin < 3
        numBootstrap = 1000; % Default number of bootstrap samples
    end

    % Fit GPD to original data
    data = sort(data,'ascend');
    pd = fitdist(data,'Logistic');
    mu = pd.mu;
    sigma = pd.sigma;
    MLE = gamcdf(Xi,mu,sigma);

    % Bootstrap resampling for confidence intervals
    CI_val = zeros(numBootstrap, length(Xi));
    for i = 1:numBootstrap
        sample = datasample(data, length(data), 'Replace', true);
        pd = fitdist(sample,'Logistic');
        mu = pd.mu;
        sigma = pd.sigma;
        CI_val(i,:)= cdf('Logistic',Xi,mu, sigma);

    end

    % Compute the confidence intervals
    CI = prctile(CI_val, [2.5, 97.5]);

end



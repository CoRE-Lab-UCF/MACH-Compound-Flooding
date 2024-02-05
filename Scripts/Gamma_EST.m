%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This function fits a gamma distribution for the selected data and CI
%   based on bootstrapping
%
%
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
function[Xi,MLE,CI]= Gamma_EST(data,Xi,numBootstrap)
    % Check input
    if nargin < 3
        numBootstrap = 1000; % Default number of bootstrap samples
    end

    % Fit GPD to original data
    data = sort(data,'ascend');
    pd = fitdist(data,'Gamma');
    mu = pd.a;
    sigma = pd.b;
    MLE = gamcdf(Xi,mu,sigma);

    % Bootstrap resampling for confidence intervals
    CI_val = zeros(numBootstrap, length(Xi));
    for i = 1:numBootstrap
        sample = datasample(data, length(data), 'Replace', true);
        pd = fitdist(sample,'Gamma');
        mu = pd.a;
        sigma = pd.b;
        CI_val(i,:)= gamcdf(Xi,mu,sigma);

    end

    % Compute the confidence intervals
    CI = prctile(CI_val, [2.5, 97.5]);

end


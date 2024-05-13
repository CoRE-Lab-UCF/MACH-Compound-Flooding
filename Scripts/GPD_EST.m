%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Pravin
%
%   This function fits a GPD distribution for the selected data and CI
%   based on bootstrapping
%
%
%   IMPORTANT:  The paths included in the script are according to the
%   author's directory. Please change them accordingly
%               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% Important: data shound be [Sample, threshold]
function[Xi,MLE,CI]= GPD_EST(data,Xi,th,numBootstrap)
    % Check input
    if nargin < 3
        numBootstrap = 1000; % Default number of bootstrap samples
    end

    % Fit GPD to original data
    paramEsts = gpfit(data,0.05);
    MLE = gpcdf(Xi,paramEsts(1),paramEsts(2),th);

    % Bootstrap resampling for confidence intervals
    CI_val = zeros(numBootstrap, size(Xi,2));
    for i = 1:numBootstrap
        sample = datasample(data, length(data), 'Replace', true);
        parm= gpfit(sample,0.05);
        CI_val(i,:) = gpcdf(Xi,parm(1),parm(2),th);
    end

    % Compute the confidence intervals
    CI = prctile(CI_val, [2.5, 97.5]);

end
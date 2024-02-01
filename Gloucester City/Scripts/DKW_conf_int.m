%This estimate the estimate and CI values for the given values in XI using
%sample data using Dvoretzky–Kiefer–Wolfowitz inequality

% Important: data shound be [Sample-threshold]
function[CI_U,CI_L,Xi]= DKW_conf_int(data,alpha)
    %sorting
    sort(data,'ascend');
    %estmiating eta
    eta = sqrt((log(2/alpha))/(2*length(data)));

    %estimating Fn(x)
    Fn= (1:1:length(data))./length(data);

    % Estimating upper bound
    CI_U = Fn+eta;
    CI_U(CI_U>1)=1;

    % Estimating lower bound
    CI_L = Fn-eta;
    CI_L(CI_L<0)=0;

    Xi=sort(data,'ascend');
end
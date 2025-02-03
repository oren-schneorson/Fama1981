function [kappa_1] = kalman1(y, kappa_0)
%KALMAN1 filter innovations given y
%   Estimating the noise a in:
%   y = a + diff(eta);  % observed (the residual from regression, Fama & ..., 1979)

y = y * sqrt(3 / var(y)); % normalize y
sigma_eta_0 = 1;
sigma_eta = 1;

sigma_a_0 = (sigma_eta_0 / kappa_0)^2 - 2*sigma_eta_0;
Var_err_y = (1+kappa_0^2)*(sigma_a_0^2+2*sigma_eta_0^2)+kappa_0*sigma_eta_0^2;


% conditionally normal
E_t_eta = @(y_) kappa_0*y_;

y_hat = NaN(T,1);
eta_hat = NaN(T+1,1);
eta_hat(1) = 0; % assume eta_0 = 0;

for t = 1:T
    y_hat(t) = -eta_hat(t); % expectation of y_t given info(t-1)
    eta_hat(t+1) = E_t_eta(y(t));
end

% compute errors
err_y = y-y_hat; % posterior errors

Var_err_y_hat = var(err_y); % variance of posterior errors
Var_err_y % theoretical variance given prior
Var_err_y_hat/Var_err_y

% get relation between real kappa and observed: Var_err_y_hat/Var_err_y

% correct kappa_0 --> kappa_1
kappa_1 = 1;





end

% TODO: finish up kalman1
% TODO: write this as a function, then go back to basic script


%{
**************************************************
    PREABMLE
**************************************************
%}

clear; clc
root = mfilename('fullpath');
username = getenv("USER")

if contains(root, 'LiveEditor') || numel(root) == 0
    % in case you're running via LiveEditor, default to root
    root = fullfile('/home', username, 'Documents', 'MATLAB', 'Fama1981');
    
else
    root = fileparts(root);
end

cd(root)

matlab_dir = fullfile('/home', username, 'Documents', 'MATLAB');
lib_data = fullfile('/media', username, 'D', 'data');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')


%rng(1)

T = 1e4;
%sigma_eta = 1.04;
%sigma_a = 1; % noise
%kappa = sigma_eta/sqrt(sigma_a+2*sigma_eta);

%kappa = 1/sqrt(3);

% kappa that would correspond to match Var_err_y == Var_err_y_hat
%kappa = .4848715643;
%kappa = .5;
precision = 1e-3;
grid_kappa = 0:precision:1;
%grid_kappa = .5;
grid_var_y = 5;
%grid_sigma_eta = 1;
%  Truth = (Uniform prior)
clf
hold on

%var_y = 3;
for var_y = grid_var_y
    %sigma_eta = sigma_eta*var(
% CANCEL THIS LOOP

ratio = NaN(size(grid_kappa));

for kappa = grid_kappa
    i = grid_kappa == kappa;

% the fact that in order to match 



sigma_eta = kappa * var_y;
% sqrt(sigma_a^2+2*sigma_eta^2) = sigma_eta / kappa
% sigma_a^2+2*sigma_eta^2 = (sigma_eta / kappa)^2
sigma_a = sqrt((sigma_eta / kappa)^2 - 2*sigma_eta^2);
%((sigma_eta / kappa)^2 - 2*sigma_eta^2 + 2*sigma_eta^2 = var_y;

a = sigma_a*randn(T, 1);
eta = sigma_eta*randn(T+1, 1);
eta(1) = 0;

y = a + diff(eta);  % observed (the residual from regression)
%y = y/var(y);


%var(y) % observed

%{
est_sigma_a = @(sigma_eta_) sqrt(var(y) - 2*sigma_eta_.^2);
est_sigma_eta = @(sigma_a_) sqrt(var(y) - sigma_a_.^2);
%}


% initial conditions
% assume sigma_a = sigma_eta; (Uniform prior)

sigma_eta_0 = 1;

%sigma_a_0 = 1;
%kappa_0 = sigma_eta_0/sqrt(sigma_a_0^2+2*sigma_eta_0^2);

%kappa_0 =.5;
%kappa_0 = .4848715643;
grid_kappa_0 = kappa * linspace(0, 5, 100);
ratio = NaN(size(grid_kappa_0));
abs_err = NaN(size(grid_kappa_0));

for kappa_0 = grid_kappa_0
    j = grid_kappa_0 == kappa_0;

%kappa_0 = kappa*2;
sigma_a_0 = (sigma_eta_0 / kappa_0)^2 - 2*sigma_eta_0;


% conditionally normal
E_t_eta = @(y_) kappa_0*y_;
V_t_eta = (1-kappa_0^2)*sigma_eta_0^2;


% expected value of y_{t+1} given y_t -- > -eta_t_hat(y_t)
E_y_t = @(y_) -E_t_eta(y_);

%var(y_{t+1} - E_y_t)


% var_eta_hat = Var(  kappa*y(t) ]  )





y_hat = NaN(T, 1);
eta_hat = NaN(T+1, 1);
eta_hat(1) = 0; % assume eta_0 = 0;

for t = 1:T
    y_hat(t) = -eta_hat(t); % expectation of y_t given info(t-1)
    eta_hat(t+1) = E_t_eta(y(t));
end

% compute errors
err_y = y-y_hat;

Var_err_y = (1+kappa_0^2)*(sigma_a_0^2+2*sigma_eta_0^2)+kappa_0*sigma_eta_0^2;
Var_err_y_hat = var(err_y);

ratio(i) = Var_err_y_hat/Var_err_y;
abs_err(i) = Var_err_y_hat;
%ratio(j) = Var_err_y_hat/Var_err_y;
%abs_err(j) = Var_err_y_hat;

end

%{
plot(grid_kappa_0, ratio)
yyaxis right
plot(grid_kappa_0, abs_err)

vline(grid_kappa(i))
vline(grid_kappa_0(ratio==max(ratio)), 'k--')
grid_kappa_0(abs_err==max(abs_err))
%}

end
plot(grid_kappa, ratio)

end


title('ratio of error variance, posterior to prior',...
    '(colors are Var(y))')
xlabel('Real kappa')
ylabel(str4fig('Var_err_y_hat / Var_err_y'))


ax = gca;

%n = numel(grid_var_y);
n = numel(grid_kappa);

cmap = winter(n);
cmap = mat2cell(cmap, ones(n, 1), 3);


[ax.Children.Color] = cmap{:};


ax = gca;

n = numel(grid_var_y);
cmap = winter(n);
colormap(cmap)
cmap = mat2cell(cmap, ones(n, 1), 3);

[ax.Children.Color] = cmap{:};

clc
colorbar

fig = gcf;
%fig.Children(1).Limits = [0 1]
%fig.Children(1).Limits = [min(grid_sigma_eta), max(grid_sigma_eta)]
%fig.Children(1).Ticks = grid_sigma_eta
fig.Children(1).Ticks = linspace(0, 1, numel(grid_var_y));
%fig.Children(1)
fig.Children(1).TickLabels = arrayfun(@num2str, grid_var_y, 'UniformOutput', false);





return


clf
subplot(1,2,1)
hold on
plot(y)
plot(y_hat)

subplot(1,2,2)
hold on
plot(eta)
plot(eta_hat)











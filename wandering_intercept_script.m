
% Wandering intercept
%{

(1) z_t = u_t - u_{t-1} + v_t  # observed signal (residual)

u_t ~ N(0, var_u)
v_t ~ N(0, var_v)

cov(u_t, u_s) = 0 for all s ~= t
cov(v_t, v_s) = 0 for all s ~= t
cov(u_t, v_s) = 0 for any s,t

get conditional expectation of u given observed z
% extract v

11-06-2024
**********************************************************
I've made significant changes to the algorithm. There was an error with
the timing of -ER, I've used -ER_t instead of -ER_{t-1}. It is fixed in
the code.

02-01-2024, Issue with sample with max date before 03/1986:
***********************************************************

The algorithm fails due to the inversion of Sigma_z, the covariance matrix
of the signal. On closer inspection, the eigenvalue of Sigma_z post
03/1986 are all above 0, but this is not the case for 02/1986 and dates
before that. possibly var_v < 0 does this...

The model assunmption is the structure in equation (1) above, along with
the iid assumptions. If those are not the case, it could cause the
calculation of var_v var_u and var_z to be negative, causing some
eigenvalue of Sigma_z to be below 0.

The autocorrelation is very high for the early period, then it goes above
-0.50 on 03/1986. This is true for the estimation of FGLS in the full
sample, but I also estimated FGLS on a limited sample (not incl 1986
onwards etc), and the result is the same - a lower than -0.50
autocorrelation of the residuals.

One solution is take the autocorrelation as a given from sample, like
1953-1999, autocorrelation of the residuals is ...
and then apply this to a shorter sample, e.g. 1953-1977.
This assumes we know the unconditional autocorrelation, and that the signal
is stationary. This might not be the case. It is likely that the standard
deviation of unexpected inflation is not constant... it is high up until
1982 and then much lower (great moderation).

I need to think how to deal with it. I'd need to have varying variance for
unexpected inflation.


%}

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


% load data, residuals from monthly estimation , FGLS
data = readtable(fullfile('/home', username, 'R', 'Fama_Gibbons_1982', 'Data', 'residulas_FGLS.xlsx'));

% going back earlier than March 1986 makes this estimation go crazy.

%min_date = datetime(1953,1,31);
%max_date = datetime(1986,3,31); 
min_date = datetime(1954,1,31);
max_date = datetime(2019,12,31);
idx_sample = data.date >= min_date & data.date <= max_date;
z = data.err_FGLS(idx_sample);
T = size(z, 1);


% variances
% ***********
months_ma = 120;
%var_z = var(z); % sample variance of signal, unconditional
var_z = arrayfun(@(t) var(z(t:t+months_ma)), 1:T-months_ma)'; % sample variance of signal, conditional
var_z = [var_z(1)*ones(months_ma, 1); var_z];
%autocorr_(z, 1) % sample 1st autocovariance
%autocorr_(z, 1) * var_z % sample estimate of cov(z_t, z_t-1)

autocorr_z = arrayfun(@(t) autocorr_(z(t:t+months_ma), 1), 1:T-months_ma)';
autocorr_z = [autocorr_z(1)*ones(months_ma, 1); autocorr_z];
%{
clf
hold on
plot(data.date(idx_sample), var_z)
yyaxis right
plot(data.date(idx_sample), autocorr_z)
legend({'var\_z', 'autocorr\_z'})

return
%}
var_u = - autocorr_(z, 1) .* var_z;
%var_u = cov(z(1:end-1), z(2:end));
%var_u = -var_u(2);

var_v = var_z - 2 * var_u;
var_v = max(var_v, 0);
var_u = (var_z-var_v)/2;



% covariances
% ***********
Sigma_u = eye(T) * var_u;

% construct Sigma_z
Sigma_z_ = 3 * eye(T) - ones(T); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
Sigma_z_ = Sigma_z_-spdiags(zeros(T, 3), -1:1, Sigma_z_);
Sigma_z = Sigma_z_ .* var_u + eye(T) .* var_v;


% construct Sigma_u_z
Sigma_u_z = 2 * eye(T) - ones(T); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
Sigma_u_z = Sigma_u_z-spdiags(zeros(T, 2), 0:1, Sigma_u_z);
Sigma_u_z = Sigma_u_z .* var_u;


%{
% if I use varying var_v and var_u, can no long invert using this
% first component, the inverse of the moving average variance-covariance matrix

% inverting Sigma_z: Sigma_z is a sum of two matrices that have closed form
% inverse matrix.
% invert a sum of knowln components inverses
% see https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices

[~, inv_Sigma_z_] = get_H_movavg(-1, T);
% second component, simple inverse of identity matrix...

% combining components following Miller, 1981 (Mathematics Magazine)
inv_Sigma_z = inv_sum(...
    Sigma_z_*var_u, eye(T)*var_v,...
    inv_Sigma_z_/var_u, eye(T)/var_v...
    );
%}

inv_Sigma_z = Sigma_z \ eye(T);



u_0 = 0;
E_u_given_z = Sigma_u_z * inv_Sigma_z * (z-0);


Delta_u = [E_u_given_z(1)-u_0; diff(E_u_given_z)];
Ev = z-Delta_u;
%a = Sigma_z\ones(T)

clf
subplot(1,2,1)
hold on
plot(data.date(idx_sample), z)
%plot(data.date, E_u_given_z)
plot(data.date(idx_sample), Ev)

subplot(1,2,2)
plot(data.date(idx_sample), -cumsum(Ev)*100)

%AR = arrayfun(@(lag) autocorr_(Ev, lag), 1:24);
%bar(AR)

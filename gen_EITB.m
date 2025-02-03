%{

In-sample procedure for TB model, FGLS estimation.

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
lib_data = fullfile('/media', username, 'D', 'data', 'BOI');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')


lib_cpi = fullfile(lib_data, 'cpi_usa');
lib_ms = fullfile(lib_data, 'money_supply_usa');
lib_ea = fullfile(lib_data, 'economic_activity');  % economic activity

lib_ir = fullfile(lib_data, 'ir');
lib_pop = fullfile(lib_data, 'population');
lib_consumption = fullfile(lib_data, 'consumption_usa');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
lag_inflation = 0;
log_approx = false;
pcent = 100; % set as 100 if want percent
freq = 'M'; % frequency to be used


%{
% original dates: Nelson 1976
sample_tab = table(...
    [1:4]',...
    [datetime(1953,6,30), datetime(1953,6,30), datetime(1953,6,30), datetime(1964,1,31)]',...
    [datetime(1971,4,30), datetime(1974,2,28), datetime(1963,12,31), datetime(1974,2,28)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}


%{
% original dates: Fama 1981
sample_tab = table(...
    [1:2]',...
    [datetime(1953,1,31), datetime(1954,1,31)]',...
    [datetime(1976,12,31), datetime(1976,2,28)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}
%
% extended dates
sample_tab = table(...
    [1:7]',...
    [datetime(1953,1,31), datetime(1953,1,31), datetime(2000,1,31),...
     datetime(1953,1,31), datetime(1953,1,31),...
     datetime(1954,1,31), datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(1977,12,31),...
    datetime(1976,12,31), datetime(1977,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%{
% paper dates
sample_tab = table(...
    [1:4]',...
    [datetime(1953,1,31), datetime(1953,1,31),...
     datetime(1953,1,31), datetime(2000,1,31)]',...
    [datetime(2019,12,31), datetime(2022,12,31),...
     datetime(1999,12,31), datetime(2019,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}


%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

table1b = [];


series_CPI = 'CPIAUCSL';

if strcmp(freq, 'M')
    series_Rf = 'RF_FF'; % 1 month TBill, from Kenneth French website, Monthly, not annuazlied
elseif strcmp(freq, 'Q')
    series_Rf = 'TB3MS'; % 3-Month Treasury Bill Secondary Market Rate, Discount Basis, Percent, Monthly, Not Seasonally Adjusted
    lib_ir = fullfile(lib_ir, 'TBs');
elseif strcmp(freq, 'A')
    series_Rf = 'TB3MS'; % 3-Month Treasury Bill Secondary Market Rate, Discount Basis, Percent, Monthly, Not Seasonally Adjusted
    %series_Rf = 'TB1YR'; % 1-Year Treasury Bill Secondary Market Rate, Discount Basis, Percent, Monthly, Not Seasonally Adjusted
    lib_ir = fullfile(lib_ir, 'TBs');
else
    error('freq must be M, Q or A')
end

%series_ms = 'BOGMBASE';
series_ms = {'CURRCIR', 'RESBALNS'};
series_ea = 'INDPRO';
series_C = {'PCENDG', 'PCES'};
series_PD = {'PCEPINDG', 'PCEPIS'};


% src: BI, Bank of Israel old website; FAME
data = struct(...
    'series', {...
    series_CPI,... Price index of all urban consumers
    series_Rf,...
    'BOGMBASE',...
    'CURRCIR',...
    'RESBALNS',...
    series_ea,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_ir,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ea,...
    %'',...
    },...
    'src', {...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    %'',...
    },...
    'VariableNames', cell(1,6)...
    );


% declare containers
metadata = struct;
for src = unique({data.src})
    src = char(src);
    metadata.(src) = table;
end


% load data to T_ and join it 
for i = 1:size(data, 2)
    
    lib_data_series = data(i).lib_data;
    [T_, metadata_] = load_data_files(data(i).series, lib_data_series, data(i).src);
    data(i).VariableNames = T_.Properties.VariableNames(2);
    
    % robust date, start of month
    aux_date = T_.date+1;
    if all(aux_date.Day == 1)
        %T_.date = aux_date-calmonths(1);
        T_.date = aux_date;
    end
    if strcmp(data(i).series, series_CPI)
        % Boons lag inflation one period, but Fama does not.
        T_.date = T_.date + calmonths(lag_inflation);
    end
    
    if ~exist('T', 'var')
        T = T_;
        metadata = metadata_;
    else
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true, 'Type', 'left');
        metadata = [metadata; metadata_];
    end
    

end


if strcmp(freq, 'M')
    %pass
elseif strcmp(freq, 'Q')
    T = table2timetable(T);
    T = retime(T, 'quarterly', 'lastvalue');
    T = timetable2table(T);
elseif strcmp(freq, 'A')
    T = table2timetable(T);
    T = retime(T, 'yearly', 'lastvalue');
    T = timetable2table(T);
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

T = fillmissing(T, 'previous');

head(T)
%help retime



params = containers.Map();
params('constant') = true;
params('weight') = false;
params('freq') = freq;



grid_h = 1:12;
% Sample period
for period = sample_tab.period(2:end)'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);

    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    idx = ismember({data.series}, series_Rf);
    T.Rf = T{:, data(idx).VariableNames}/100*pcent;
    if strcmp(freq, 'M')
        % pass
    elseif strcmp(freq, 'Q')
        T.Rf = T.Rf/4;
    elseif strcmp(freq, 'A')
        % pass
    else
        error('freq must be M, Q or A')
    end

    
    T.Rf_lag = lag(T.Rf, 1);
    T.Rf_lead = lead(T.Rf, 1);
    
    idx = ismember({data.series}, series_CPI);
    T.Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    
    
    if log_approx
        T.Pi = log(T.Pi)*pcent;
    else
        T.Pi = (T.Pi -1)*pcent;
    end

    T_ = T(idx_sample, :);

    summary_stats = [...
        [arrayfun(@(h) autocorr_(T_.Pi, h), grid_h), mean(T_.Pi), std(T_.Pi)];...
        [arrayfun(@(h) autocorr_(T_.Rf, h), grid_h), mean(T_.Rf), std(T_.Rf)]...
        ];
    array2table(summary_stats, 'VariableNames',...
        [strcat('rho_', arrayfun(@num2str, grid_h, 'UniformOutput', false)),...
        {'mean', 'std'}])

    Y = diff(T_.Pi);
    X = diff(T_.Rf);
    
    %X = [ones(size(X)), X];

    N = size(X, 1);
    K = size(X, 2);

    H_hat = @(theta_) get_H_movavg(theta_, N);
    Omega_inv_hat_ = @(theta_) H_hat(theta_)'*H_hat(theta_);

    beta_FGLS_fun = @(theta_) ((X'*Omega_inv_hat_(theta_)*X)\eye(K))*...
        X'*Omega_inv_hat_(theta_)*Y;
    sigma2_fun = @(beta_) (Y-X*beta_)'*(Y-X*beta_)/N;
    Omega_inv_hat = @(theta_) Omega_inv_hat_(theta_)/sigma2_fun(beta_FGLS_fun(theta_));

    %beta_FGLS = @(theta) X'*Omega_inv_hat(theta)*Y;
    % observe that errors are weighted by H_hat

    err = @(theta_) (Y-X*beta_FGLS_fun(theta_));
    errw = @(theta_) H_hat(theta)*(Y-X*beta_FGLS_fun(theta_));
    J = @(theta_) sqrt(sum(err(theta_).^2, 'omitnan')/N);

    % don't minimize SSR for theta, that is inconsistent. 
    % use *log*-likelihood
    % implement maximum likelihood with and without FGLS assumption.

    % log-likelihood function
    L_fun = @(theta_, beta_)...
        -N/2*log(2*pi)-.5*sum(log(...
        eig(get_Omega_movavg(theta_, sigma2_fun(beta_), N))))-...
        (Y-X*beta_)'*Omega_inv_hat(theta_)*(Y-X*beta_)/2/sigma2_fun(beta_);

    % log-likelihood function, assuming beta=beta_FGLS
    L_fun_FGLS = @(theta_)...
        -N/2*log(2*pi)-.5*sum(log(...
        eig(get_Omega_movavg(theta_, sigma2_fun(beta_FGLS_fun(theta_)), N))))-...
        (Y-X*beta_FGLS_fun(theta_))'*...
        Omega_inv_hat(theta_)*...
        (Y-X*beta_FGLS_fun(theta_))/2;

    %%
    % Maximum likeloihhod estimation, jointly on (theta, beta); iteratively
    % *******************************************************************
    clc
    % init conditions
    theta_prev = -1;
    theta_ML = -.8;
    precision = 1e-5;
    precision_theta = 1e-2;
    precision_beta = 1e-2;
    spectrum_theta = [-1, -.2];
    spectrum_beta = [.1, 2];
    while abs(theta_ML-theta_prev) > precision
        theta_prev = theta_ML;
        
        [beta_ML, ~] = my_fminunc(@(beta_) -L_fun(theta_ML, beta_), spectrum_beta, precision_beta);
        [theta_ML, ~] = my_fminunc(@(theta_) -L_fun(theta_, beta_ML), spectrum_theta, precision_theta);
        
        precision_theta = max(precision_theta/10, 1e-5);
        precision_beta = max(precision_beta/10, 1e-5);

        spectrum_beta = [beta_ML-10*precision_beta, beta_ML+10*precision_beta];
        spectrum_theta = [theta_ML-10*precision_theta, theta_ML+10*precision_theta];
        spectrum_theta = max(spectrum_theta, -1);

        fprintf('Fixed point estimation: theta_{ML}=%.4f  beta_{ML}=%.4f\n',...
            theta_ML, beta_ML)
    end

    % Maximum likeloihhod estimation, on theta + FGLS assumption; iteratively
    % *******************************************************************
    clc
    % init conditions
    theta_prev = -1;
    theta_FGLS = -.8;
    precision = 1e-5;
    precision_theta = 1e-2;
    spectrum_theta = [-1, -.2];

    while abs(theta_FGLS-theta_prev) > precision
        theta_prev = theta_FGLS;
        
        [theta_FGLS, ~] = my_fminunc(@(theta_) -L_fun_FGLS(theta_), spectrum_theta, precision_theta);
        beta_FGLS = beta_FGLS_fun(theta_FGLS);
        
        precision_theta = max(precision_theta/10, 1e-5);
        spectrum_theta = [theta_FGLS-10*precision_theta, theta_FGLS+10*precision_theta];
        spectrum_theta = max(spectrum_theta, -1);
        
        fprintf('Fixed point estimation: \theta_{FGLS}=%.4f  \beta_{FGLS}=%.4f\n',...
            theta_FGLS, beta_FGLS)

    end
    
    %{
    %% plot the joint log-likelihood function

    grid_theta = linspace(-1, -.5, 26);
    grid_beta = linspace(.5, 1.5, 25);

    grid_theta2 = repmat(grid_theta, numel(grid_beta), 1)'; % ax=1
    grid_beta2 = repmat(grid_beta, numel(grid_theta), 1); % ax=2
    l = arrayfun(@(i) L_fun(grid_theta2(i), grid_beta2(i)), 1:numel(grid_theta2));
    l = reshape(l, [numel(grid_theta), numel(grid_beta)]);
    
    % check reshape
    %all(all(reshape(grid_theta2(:), [numel(grid_theta), numel(grid_beta)])==grid_theta2))
    %all(all(reshape(grid_beta2(:), [numel(grid_theta), numel(grid_beta)])==grid_beta2))
    
    clf
    hold on
    contourf(grid_theta, grid_beta, l')
    colorbar
    scatter(theta_ML, beta_ML, 'r*')
    scatter(theta_FGLS, beta_FGLS, 'k*')
    xlabel('\theta (moving average parameter)')
    ylabel('\beta')
    title('Joint log-likelihood function',...
        sprintf('United States, %s -- %s', min_date, max_date))
    legend({'', 'Maximum likelihood', 'FGLS'})
    return
    %%
    %}
    
    [beta_OLS, ~, ~, ~, STATS] = regress(Y, X);

    Y_hat_ML = X*beta_ML;
    Y_hat_FGLS = X*beta_FGLS;
    Y_hat_OLS = X*beta_OLS;
    Y_hat_FAMA = X * 1.15; % only relevant for 1953-1977
    Y_hat_1 = X * 1.00;

    err_ML = Y-Y_hat_ML;
    err_FGLS = Y-Y_hat_FGLS;
    err_OLS = Y-Y_hat_OLS;
    err_FAMA = Y-Y_hat_FAMA;
    err_1 = Y-Y_hat_1;
    
    errw_ML = H_hat(theta_ML)*(Y-Y_hat_ML);
    errw_FGLS = H_hat(theta_FGLS)*(Y-Y_hat_FGLS);
    errw_FAMA = H_hat(theta_FGLS)*(Y-Y_hat_FAMA);

    SSR_ML_diff = sum(err_ML.^2);
    SSR_FGLS_diff = sum(err_FGLS.^2);
    SSR_FAMA_diff = sum(err_FAMA.^2);
    SSR_1_diff = sum(err_1.^2);
    SS_tot_diff = sum((Y-mean(Y)).^2);%for period = sample_tab.period(5)'


    R_sq_ML_diff = 1-SSR_ML_diff/SS_tot_diff;
    R_sq_FGLS_diff = 1-SSR_FGLS_diff/SS_tot_diff;
    R_sq_OLS_diff = STATS(1);
    R_sq_FAMA_diff = 1-SSR_FAMA_diff/SS_tot_diff;
    R_sq_1_diff = 1-SSR_1_diff/SS_tot_diff;

    %{
    %%
    clf
    clc
    [min_date, max_date]
    [beta_ML, beta_FGLS, beta_OLS, 1.15, 1]
    [R_sq_ML_diff, R_sq_FGLS_diff, R_sq_OLS_diff, R_sq_FAMA_diff, R_sq_1_diff]

    hold on
    plot(T_.date(2:end), err_FGLS)
    plot(T_.date(2:end), err_1)
    %{
    1953-1977, same result as FAMA, R-sq is not very good in the difference
    equation. Perhaps this is to be expected, because the standard
    deviation of diff(Rf) is much lower than standard deviation of
    inflation.
    %}
    
    return
    %%
    %}

    % data to save for R to run a KalmanFilter procedure
    fname_R_in = sprintf('%s--%s_residulas_ML_FGLS.xlsx', min_date, max_date);
    fpath_R_in = fullfile('/home', username, 'R', 'Fama_Gibbons_1982', 'data', 'us', fname_R_in);
    writetable([T_(2:end, 'date'), array2table(err_ML), array2table(err_FGLS)],...
        fpath_R_in)
    
    % output of KalmanFilter procedure from R
    fname_R_out = sprintf('%s--%s_KalmanFilter.csv', min_date, max_date);
    fpath_R_out = fullfile('/home', username, 'R', 'Fama_Gibbons_1982', 'results', 'us', fname_R_out);
    R_out = readtable(fpath_R_out);
    R_out.Properties.VariableNames = {'date', 'u_t', 'u_t_', 'v_t'};
    R_out.date = T_.date;
    %%
    %{
    ***************************************************
    Wandering intercept procedure
    ***************************************************
    code extract from wandering_intercept.m, uses residulas_FGLS.xlsx
    

    % 
    z_t = u_t - u_{t-1} + v_t  # observed signal (residual)
    
    u_t ~ N(0, var_u)
    v_t ~ N(0, var_v)
    
    cov(u_t, u_s) = 0 for all s ~= t
    cov(v_t, v_s) = 0 for all s ~= t
    cov(u_t, v_s) = 0 for any s,t
    
    get conditional expectation of u given observed z
    % extract v (conditional expectation...)
    
    %}
        
    z_ML = err_ML; % set the observed signal as err_FGLS
    z_FGLS = err_FGLS; % set the observed signal as err_FGLS
    z_1 = err_1; % set the observed signal as err_FGLS
    
    % variances
    % ***********
    % note: using unconditional var_z, var_v and var_u causes
    % complications and is economically unsound. See wandering_intercept.m
    years_ma = 5;
    if strcmp(freq, 'M')
        periods_ma = years_ma*12;
    elseif strcmp(freq, 'Q')
        periods_ma = years_ma*4;
    elseif strcmp(freq, 'A')
        periods_ma = years_ma;
    else
        error('freq must be M, Q or A')
    end

    var_z_ML = var(z_ML); % sample variance of signal, unconditional
    %var_z_ML = arrayfun(@(t) var(z_ML(t:t+periods_ma)), 1:N-periods_ma)'; % sample variance of signal, conditional
    %var_z_ML = [var_z_ML(1)*ones(periods_ma, 1); var_z_ML];

    
    var_z_1 = var(z_1); % sample variance of signal, unconditional (beta=1)
    %var_z_1 = arrayfun(@(t) var(z_1(t:t+periods_ma)), 1:N-periods_ma)'; % sample variance of signal, conditional
    %var_z_1 = [var_z_1(1)*ones(periods_ma, 1); var_z_1];

    % first order autocorrelation
    autocorr_z_ML = autocorr_(z_ML, 1);
    %autocorr_z_ML = arrayfun(@(t) autocorr_(z_ML(t:t+periods_ma), 1), 1:N-periods_ma)';
    %autocorr_z_ML = [autocorr_z_ML(1)*ones(periods_ma, 1); autocorr_z_ML];

    autocorr_z_1 = autocorr_(z_1, 1);
    %autocorr_z_1 = arrayfun(@(t) autocorr_(z1(t:t+periods_ma), 1), 1:N-periods_ma)';
    %autocorr_z_1 = [autocorr_z_1(1)*ones(periods_ma, 1); autocorr_z_1];

    % cov(z_t, z_{t-1}) = -var_u
    % autocorr_z = cov(z_t, z_{t-1}) / var_z
    var_u_ML = - autocorr_z_ML .* var_z_ML;
    var_u_1 = - autocorr_z_1 .* var_z_1;

    var_v_ML = var_z_ML - 2 * var_u_ML;
    var_v_1 = var_z_1 - 2 * var_u_1;

    if var_v_ML < 0
        warning('Sample autocorrelation of signal lower than -0.50. Setting var_v to 10% of var_z')
        var_v_ML = 0.1*var_z_ML;
        var_v_1 = 0.1*var_z_1;
        var_u_ML = (var_z_ML-var_v_ML)/2;
        var_u_1 = (var_z_1-var_v_1)/2;

    elseif var_u_ML < 0
        warning('Sample autocorrelation of signal greater than 0.00. Setting var_v to 90% of var_z')
        var_v_ML = 0.9*var_z_ML;
        var_v_1 = 0.9*var_z_1;
        var_u_ML = (var_z_ML-var_v_ML)/2;
        var_u_1 = (var_z_1-var_v_1)/2;

    end

    %var_v_ML = max(var_v_ML, 0);
    %var_v_1 = max(var_v_1, 0);
    %var_u_ML = (var_z_ML-var_v_ML)/2;
    %var_u_1 = (var_z_1-var_v_1)/2;

    %autocorr_(z_ML, 1) % sample 1st autocovariance
    %autocorr_(z_ML, 1) * var_z % sample estimate of cov(z_t, z_t-1)
    %var_z_ML/var_v_ML % signal to noise ratio
    
    
    % covariances
    % ***********
    Sigma_u = eye(N) .* var_u_ML;
    Sigma_u1 = eye(N) .* var_u_1;
    
    % construct Sigma_z
    Sigma_z_ = 3 * eye(N) - ones(N); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
    Sigma_z_ = Sigma_z_-spdiags(zeros(N, 3), -1:1, Sigma_z_);
    Sigma_z = Sigma_z_ .* var_u_ML + eye(N) .* var_v_ML;

    % construct Sigma_z1
    Sigma_z1_ = 3 * eye(N) - ones(N); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
    Sigma_z1_ = Sigma_z1_-spdiags(zeros(N, 3), -1:1, Sigma_z1_);
    Sigma_z1 = Sigma_z1_ .* var_u_1 + eye(N) .* var_v_1;

    % construct Sigma_u_z
    Sigma_u_z = 2 * eye(N) - ones(N); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
    Sigma_u_z = Sigma_u_z-spdiags(zeros(N, 2), 0:1, Sigma_u_z);
    Sigma_u_z = Sigma_u_z .* var_u_ML;

    % construct Sigma_u_z1
    Sigma_u_z1 = 2 * eye(N) - ones(N); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
    Sigma_u_z1 = Sigma_u_z1-spdiags(zeros(N, 2), 0:1, Sigma_u_z1);
    Sigma_u_z1 = Sigma_u_z1 .* var_u_1;


    if numel(var_z_ML) ~= 1
        % if I use varying var_v and var_u, can no longer invert using this
        % first component, the inverse of the moving average
        % variance-covariance matrix

        inv_Sigma_z = Sigma_z \ eye(N);
        inv_Sigma_z1 = Sigma_z1 \ eye(N);

    else
        % inverting Sigma_z: Sigma_z is a sum of two matrices that have closed form
        % inverse matrix.
        % invert a sum of knowln components inverses
        % see https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
        
        % first component, the inverse of the moving average variance-covariance matrix
        [~, inv_Sigma_z_] = get_H_movavg(-1, N);
        % second component, simple inverse of identity matrix...
        
        % combining components following Miller, 1981 (Mathematics Magazine)
        inv_Sigma_z = inv_sum(...
            Sigma_z_*var_u_ML, eye(N)*var_v_ML,...
            inv_Sigma_z_/var_u_ML, eye(N)/var_v_ML...
            );
        inv_Sigma_z1 = inv_sum(...
            Sigma_z1_*var_u_1, eye(N)*var_v_1,...
            inv_Sigma_z_/var_u_1, eye(N)/var_v_1...
            );
    end
    
    u_0 = 0;
    
    % note: E_u|z has a lot of autocorrelation, which translates later to
    % UITB. This should not happen, UITB is to be clear of autocorrelation.
    % Is this because it is expectation? Is it solved by using a kalman
    % filter? Or rather, I should use this only OOS?
    E_u_given_z_ML = Sigma_u_z * inv_Sigma_z * (z_ML-0);
    Delta_u_ML = [E_u_given_z_ML(1)-u_0; diff(E_u_given_z_ML)];
    Ev_ML = z_ML-Delta_u_ML;

    E_u_given_z1 = Sigma_u_z1 * inv_Sigma_z1 * (z_1-0);
    Delta_u_1 = [E_u_given_z1(1)-u_0; diff(E_u_given_z1)];
    Ev_1 = z_1-Delta_u_1;

    %
    % Comparison of KalmanFilter method with CondNorm method
    % ******************************************************
    %%
    clc
    acf_err_ML = arrayfun(@(h) autocorr_(err_ML, h), grid_h);
    acf_err_FGLS = arrayfun(@(h) autocorr_(err_FGLS, h), grid_h);
    acf_err_OLS = arrayfun(@(h) autocorr_(err_OLS, h), grid_h);
    acf_errw_FGLS = arrayfun(@(h) autocorr_(errw_FGLS, h), grid_h);
    acf_E_u_given_z = arrayfun(@(h) autocorr_(E_u_given_z_ML, h), grid_h);
    acf_R_out_u_t = arrayfun(@(h) autocorr_(R_out.u_t(2:end), h), grid_h);

    
    clf
    subplot(1,2,1)
    hold on
    plot(R_out.date(2:end), R_out.u_t(2:end), 'LineWidth', 2)
    plot(T_.date(2:end), E_u_given_z_ML, '--')

    corr_ = corrcoef(R_out.u_t(2:end), E_u_given_z_ML);
    corr_ = corr_(2);
    title(sprintf('mean: KalmanFilter=%.3f; CondNormal=%.3f',...
        mean(R_out.u_t(2:end)), mean(E_u_given_z_ML)),...
        sprintf('correlation: %.3f', corr_))
    legend({'u_t (KalmanFilter)', 'u_t (CondNormal)'})

    subplot(1,2,2)
    %bar(grid_h, [acf_u_FGLS; acf_u_FAMA; acf_err_OLS; acf_err_FGLS])
    %legend({'u_{FGLS}', 'u_{FAMA}', 'err_{OLS}', 'err_{FGLS}'})
    bar(grid_h, [acf_E_u_given_z; acf_R_out_u_t])
    title(sprintf('ML: ACF(1)'))
    legend({'u_t (KalmanFilter)', 'u_t (CondNormal)'})
    ylim([-1 1])
    fname_comparing_KalmanFilter2CondNormal = sprintf('US__%s--%s.jpg',...
        min_date, max_date);
    fpath_comparing_KalmanFilter2CondNormal = fullfile('.', 'gfx',...
        'comparing_KalmanFilter2CondNormal', 'us',...
        fname_comparing_KalmanFilter2CondNormal);
    export_fig(fpath_comparing_KalmanFilter2CondNormal)
    continue
    %%
    %}

    %{
    %%
    clc
    clf
    plot(T_.date(2:end), -cumsum(Ev)*100)
    %%
    %}
    resid = T_.Pi-beta_ML*T_.Rf - [0; cumsum(Ev_ML)];
    alpha = [0; cumsum(Ev_ML)] + mean(resid);

    resid1 = T_.Pi-1*T_.Rf - [0; cumsum(Ev_1)];
    alpha1 = [0; cumsum(Ev_1)] + mean(resid1);

    EITB = alpha(1:end-1) + beta_FGLS * T_.Rf(2:end);
    UITB = T_.Pi(2:end)-EITB;
    EITB1 = alpha1(1:end-1) + 1.0 * T_.Rf(2:end);
    UITB1 = T_.Pi(2:end)-EITB1;
    
    clf
    subplot(1,2,1)
    hold on
    plot(T_.date(2:end), T_.Pi(2:end)*100/pcent)
    plot(T_.date(2:end), EITB*100/pcent)
    %plot(T_.date(2:end), EITB1*100/pcent, '--')
    subplot(1,2,2)
    hold on
    bar(arrayfun(@(h) autocorr_(UITB, h), 1:24))
    %bar(arrayfun(@(h) autocorr_(UITB1, h), 1:24))


    %plot(T_.date(2:end), UITB)
    %plot(T_.date(2:end), UITB1)
    % plot moving average
    clf
    hold on
    h=0;
    plot(T_.date(2+h:end), arrayfun(@(t) mean(UITB(t:t+h)), 1:N-h)')
    h=12;
    plot(T_.date(2+h:end), arrayfun(@(t) mean(UITB(t:t+h)), 1:N-h)')

    drawnow

    SSR = sum(UITB.^2);
    SSR_1 = sum(UITB1.^2);
    SS_tot = sum((T_.Pi(3:end)-mean(T_.Pi(3:end))).^2);
    R_sq = 1-SSR/SS_tot
    R_sq_1 = 1-SSR_1/SS_tot
    %}

    [min_date, max_date]
    [1, beta_FGLS]
    [R_sq_1, R_sq]
    
    return
    var(UITB)/var(z_ML)
    (var(z_ML)-2*var(UITB))/var(z_ML)

    EITB = [T_(2:end, 'date'),...
        array2table(EITB*100/pcent, 'VariableNames', {'EITB'})];
    UITB = [T_(2:end, 'date'),...
        array2table(UITB*100/pcent, 'VariableNames', {'UITB'})];
    EITB1 = [T_(2:end, 'date'),...
        array2table(EITB1*100/pcent, 'VariableNames', {'EITB1'})];
    UITB1 = [T_(2:end, 'date'),...
        array2table(UITB1*100/pcent, 'VariableNames', {'UITB1'})];

    
    fname = sprintf('EITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB = fullfile(root, 'EITB', fname);
    writetable(EITB, fpath_EITB)
    
    fname = sprintf('UITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB = fullfile(root, 'UITB', fname);
    writetable(UITB, fpath_UITB)

    fname = sprintf('EITB1__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB1 = fullfile(root, 'EITB', fname);
    writetable(EITB1, fpath_EITB1)
    
    fname = sprintf('UITB1__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB1 = fullfile(root, 'UITB', fname);
    writetable(UITB1, fpath_UITB1)

    
end



return

%%

period = 1;
min_date = sample_tab.min_date(period);
max_date = sample_tab.max_date(period);


fpath = fullfile(matlab_dir, 'NRC', 'arima_101_usa.csv')
arima_101 = readtable(fpath);
arima_101 = renamevars(arima_101, 'Pi_hat', 'Pi_hat_101');
arima_101.date = arima_101.date+1;
arima_101.date = arima_101.date - calmonths(1);
arima_101.date = arima_101.date-1;

fname = sprintf('EITB__%s-%s.%s.csv', min_date, max_date, freq);
fpath_EITB = fullfile(root, 'EITB', fname);
fname = sprintf('UITB__%s-%s.%s.csv', min_date, max_date, freq);
fpath_UITB = fullfile(root, 'UITB', fname);

EITB = readtable(fpath_EITB);
UITB = readtable(fpath_UITB);

TB_model = join(EITB, UITB);
TB_model
arima_101

common_dates = intersect(arima_101.date, TB_model.date);
idx_arima = ismember(arima_101.date, common_dates);
idx_TB = ismember(TB_model.date, common_dates);
idx_T = ismember(T.date, common_dates);

clf
hold on
%plot(T.date(idx_T), T.Pi(idx_T), 'r-')
plot(arima_101.date(idx_arima), arima_101.Pi_hat_101(idx_arima), 'k-')
plot(TB_model.date(idx_TB), TB_model.EITB(idx_TB), 'r--')

sqrt(mean((arima_101.epsilon_101(idx_arima)).^2))
sqrt(mean((TB_model.UITB(idx_TB)).^2))

%legend({'\pi', 'EIAR', 'EITB'})
legend({'EIAR', 'EITB'})
ylabel('Percent')


corr(TB_model.EITB(idx_TB), arima_101.Pi_hat_101(idx_arima))
corr(TB_model.UITB(idx_TB), arima_101.epsilon_101(idx_arima))





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


matlab_dir = fullfile('/home', username, 'Documents', 'MATLAB');
lib_data = fullfile('/media', username, 'D', 'data');

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
lag_inflation = 1;
log_approx = false;
pcent = 1; % set as 100 if want percent
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
for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);



    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    idx = ismember({data.series}, series_Rf);
    T.Rf = T{:, data(idx).VariableNames}/100;
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

    Y = diff(T_.Pi(2:end));
    X = diff(T_.Rf(2:end));
    %X = [ones(size(X)), X];

    N = size(X, 1);
    K = size(X, 2);

    H_hat = @(theta) get_H_movavg(theta, N);
    Omega_inv_hat = @(theta) H_hat(theta)'*H_hat(theta);
    beta_FGLS_fun = @(theta) ((X'*Omega_inv_hat(theta)*X)\eye(K))*...
        X'*Omega_inv_hat(theta)*Y;
    %beta_FGLS = @(theta) X'*Omega_inv_hat(theta)*Y;
    % observe that errors are weighted by H_hat
    err = @(theta) H_hat(theta)*(Y-X*beta_FGLS_fun(theta));
    J = @(theta) sum(err(theta).^2, 'omitnan')/N;


    step = 1e-2; % precision parameter for theta_estimation
    theta_grid = -1:step:-.5;
    res = arrayfun(J, theta_grid);
    theta_min = theta_grid(res == min(res));

    clf
    plot(theta_grid, res)
    vline(theta_grid(res==min(res)))
    drawnow
    
    step = 1e-3; % precision parameter for theta_estimation
    theta_grid = [theta_min-step*25:step:theta_min+step*25];
    res = arrayfun(J, theta_grid);
    theta_min = theta_grid(res == min(res));

    clf
    plot(theta_grid, res)
    vline(theta_grid(res==min(res)))
    drawnow

    step = 1e-4; % precision parameter for theta_estimation
    theta_grid = [theta_min-step*25:step:theta_min+step*25];
    res = arrayfun(J, theta_grid);

    clf
    plot(theta_grid, res)
    vline(theta_grid(res==min(res)))
    drawnow

    step = 1e-5; % precision parameter for theta_estimation
    theta_grid = [theta_min-step*25:step:theta_min+step*25];
    res = arrayfun(J, theta_grid);

    clf
    plot(theta_grid, res)
    vline(theta_grid(res==min(res)))
    drawnow

    %{
    clf
    plot(theta_grid, res)
    % result is very close to FAMA
    return
    %}

    % compare estimators

    clc
    theta = theta_grid(res == min(res));
    beta_FGLS = beta_FGLS_fun(theta);
    [beta_OLS, ~, ~, ~, STATS] = regress(Y, X);

    % STATS(1) % R_sq
    Y_hat_FGLS = X*beta_FGLS_fun(theta);
    Y_hat_OLS = X*beta_OLS;
    Y_hat_FAMA = X * 1.15;
    %Y_hat_FAMA = X * 1.00;

    err_OLS = Y-Y_hat_OLS;
    err_FGLS = Y-Y_hat_FGLS;
    err_FAMA = Y-Y_hat_FAMA;

    
    errw_FGLS = H_hat(theta)*(Y-Y_hat_FGLS);
    errw_FAMA = H_hat(theta)*(Y-Y_hat_FAMA);
    grid_h = 1:12;

    %sum(err_FGLS.^2)/N
    %sum(err_FAMA.^2)/N
    

    % err_FGLS_t = u_t - theta * u_{t-1}
    u_0 = 0;
    u_FGLS = errw_FGLS + theta * [u_0; errw_FGLS(1:end-1)];
    u_FAMA = errw_FAMA + theta * [u_0; errw_FAMA(1:end-1)];


    %{
    acf_err_OLS = arrayfun(@(h) autocorr_(err_OLS, h), grid_h);
    acf_err_FGLS = arrayfun(@(h) autocorr_(err_FGLS, h), grid_h);
    acf_errw_FGLS = arrayfun(@(h) autocorr_(errw_FGLS, h), grid_h);
    acf_u_FGLS = arrayfun(@(h) autocorr_(u_FGLS, h), grid_h);
    acf_u_FAMA = arrayfun(@(h) autocorr_(u_FAMA, h), grid_h);
    clf
    %bar(grid_h, [acf_u_FGLS; acf_u_FAMA; acf_err_OLS; acf_err_FGLS])
    %legend({'u_{FGLS}', 'u_{FAMA}', 'err_{OLS}', 'err_{FGLS}'})
    bar(grid_h, [acf_err_FGLS; acf_errw_FGLS])
    legend({'err_{FGLS}', 'errw_{FGLS}'})
    % here we see that Fama's residulas, and the breakdown to u are very
    % close. The autocorrelation in u is still substantial.

    %}

    writetable([T_(3:end, 'date'), array2table(err_FGLS), array2table(u_FGLS)],...
        'residulas_FGLS.xlsx')

    %{
    ***************************************************
    Wandering intercept procedure
    ***************************************************
    code extract from wandering_intercept.m, uses residulas_FGLS.xlsx
    

    
    z_t = eta_t - eta_{t-1} + v_t  # observed signal (residual)
    
    eta_t ~ N(0, var_eta)
    v_t ~ N(0, var_v)
    
    cov(eta_t, eta_s) = 0 for all s ~= t
    cov(v_t, v_s) = 0 for all s ~= t
    cov(eta_t, v_s) = 0 for any s,t
    
    get conditional expectation of eta given observed z
    % extract v (conditional expectation...)
    
    %}
    
    
    z = err_FGLS; % set the observed signal as err_FGLS
    sample_size = size(z, 1);
    
    
    % variances
    % ***********
    % note: using unconditional var_z, var_v and var_eta causes
    % complications and is economically unsound. See wandering_intercept.m
    years_ma = 10;
    if strcmp(freq, 'M')
        months_ma = years_ma*12;
    elseif strcmp(freq, 'Q')
        months_ma = years_ma*4;
    elseif strcmp(freq, 'A')
        months_ma = years_ma;
    else
        error('freq must be M, Q or A')

    end

    %var_z = var(z); % sample variance of signal, unconditional
    var_z = arrayfun(@(t) var(z(t:t+months_ma)), 1:sample_size-months_ma)'; % sample variance of signal, conditional
    var_z = [var_z(1)*ones(months_ma, 1); var_z];

    autocorr_z = arrayfun(@(t) autocorr_(z(t:t+months_ma), 1), 1:sample_size-months_ma)';
    autocorr_z = [autocorr_z(1)*ones(months_ma, 1); autocorr_z];

    var_eta = - autocorr_(z, 1) .* var_z;
    var_v = var_z - 2 * var_eta;
    var_v = max(var_v, 0);
    var_eta = (var_z-var_v)/2;

    %autocorr_(z, 1) % sample 1st autocovariance
    %autocorr_(z, 1) * var_z % sample estimate of cov(z_t, z_t-1)
    %var_z/var_v % signal to noise ratio
    
    
    % covariances
    % ***********
    Sigma_eta = eye(sample_size) .* var_eta;
    
    % construct Sigma_z
    Sigma_z_ = 3 * eye(sample_size) - ones(sample_size); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
    Sigma_z_ = Sigma_z_-spdiags(zeros(sample_size, 3), -1:1, Sigma_z_);
    Sigma_z = Sigma_z_ .* var_eta + eye(sample_size) .* var_v;
    
    % construct Sigma_eta_z
    Sigma_eta_z = 2 * eye(sample_size) - ones(sample_size); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
    Sigma_eta_z = Sigma_eta_z-spdiags(zeros(sample_size, 2), 0:1, Sigma_eta_z);
    Sigma_eta_z = Sigma_eta_z .* var_eta;
    

    %{
    % if I use varying var_v and var_eta, can no long invert using this
    % first component, the inverse of the moving average variance-covariance matrix


    % inverting Sigma_z: Sigma_z is a sum of two matrices that have closed form
    % inverse matrix.
    % invert a sum of knowln components inverses
    % see https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
    
    % first component, the inverse of the moving average variance-covariance matrix
    [~, inv_Sigma_z_] = get_H_movavg(-1, sample_size);
    % second component, simple inverse of identity matrix...
    
    % combining components following Miller, 1981 (Mathematics Magazine)
    inv_Sigma_z = inv_sum(...
        Sigma_z_*var_eta, eye(sample_size)*var_v,...
        inv_Sigma_z_/var_eta, eye(sample_size)/var_v...
        );
    %}

    inv_Sigma_z = Sigma_z \ eye(sample_size);
    
    eta_0 = 0;
    E_eta_given_z = Sigma_eta_z * inv_Sigma_z * (z-0);
    Delta_eta = [E_eta_given_z(1)-eta_0; diff(E_eta_given_z)];
    Ev = z-Delta_eta;
    
    clc
    clf
    %plot(T_.date(3:end), -cumsum(Ev)*100)

    %plot(T_.date, T_.Pi*100)
    resid = T_.Pi(3:end)-beta_FGLS_fun(theta)*T_.Rf(3:end) - cumsum(Ev);
    alpha = cumsum(Ev) + mean(resid);
    resid_1 = T_.Pi(3:end)-1.0*T_.Rf(3:end) - cumsum(Ev);
    alpha_1 = cumsum(Ev) + mean(resid_1);

    EITB = alpha + beta_FGLS_fun(theta) * T_.Rf(3:end);
    UITB = T_.Pi(3:end)-EITB;
    EITB_1 = alpha_1 + 1.0 * T_.Rf(3:end);
    UITB_1 = T_.Pi(3:end)-EITB_1;
    
    clf
    hold on
    plot(T_.date(3:end), T_.Pi(3:end)*100)
    plot(T_.date(3:end), EITB*100)
    %plot(T_.date(3:end), UITB)
    %plot(T_.date(3:end), UITB_1)

    drawnow
    return
    SSR = sum(UITB.^2);
    SSR_1 = sum(UITB_1.^2);
    SS_tot = sum((T_.Pi(3:end)-mean(T_.Pi(3:end))).^2);
    R_sq = 1-SSR/SS_tot
    R_sq_1 = 1-SSR_1/SS_tot
    %}

    EITB = [T_(3:end, 'date'),...
        array2table(EITB*100, 'VariableNames', {'EITB'})];
    UITB = [T_(3:end, 'date'),...
        array2table(UITB*100, 'VariableNames', {'UITB'})];
    
    
    fname = sprintf('EITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB = fullfile(root, 'EITB', fname);
    writetable(EITB, fpath_EITB)
    
    fname = sprintf('UITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB = fullfile(root, 'UITB', fname);
    writetable(UITB, fpath_UITB)

    
end




%%

period = 1;
min_date = sample_tab.min_date(period);
max_date = sample_tab.max_date(period);


fpath = fullfile(matlab_dir, 'NRC', 'arima_101_usa.csv');
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





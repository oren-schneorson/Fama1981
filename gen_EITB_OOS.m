%{

This script generates FGLS estimates and decomposition of autocorrelated
residual into moving average + white noise -- for each t, using information
up to t.

the model is too good. I need to withold date t observation, and attempt a
real OOS performance.

Need to think hard about timing of Rf and Pi.

Pi:
The data first reads inflation as the inflation actualized in date t.
I then lag inflation one month, to account for the delay in information.
For example, if inflation was 0.1% in between 1-31 Oct 2023, this is observed only 
on 12 Dec 2023, a lag of about 6 weeks.The observation is dated in FRED
01/11/2023

Interest rate data is attained from FRED, which has start of month dates. I
checked, and it corresponds to end of last month, so minus one day.

Thus, when I attempt to predict inflation on Nov 2023, I need the interest
rate that was in effect on 30 Nov 2023 (two weeks, most recent interest
rate before information revealed. Then TB is at t-1 to inflation.
In effect, by lagging inflation one month, I achieve this already. Let's
see why.

Inflation between 1-30 Oct 2023 is dated 01/11/2023, 31/10/2023 here.
I lag this one month, so it is dated 30/11/2023 (t). TBill rates are dated
30/11/2023 (t). I run a model using data up to t-1, i.e. both interest rates
and TBill rates that are available up until 31/10/2023 (t-1). This gives me the
parameters of the model (incl beta_FGLS and the wandering intercept).
Finally, I use the observation of 30/11/2023 (t) of the TBill rate + the
parameters of the model estimated using information up to t-1, to predict
the time t observation of inflation (dated 30/11/2023 in my data, in fact
published on 12/12/2023...



Now, for the estimation, I need to evaluate the model up to t-1, and then
guess inflation at t.

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
% I'm using first of month here, then later minus one day.
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


params = containers.Map();
params('constant') = true;
params('weight') = false;
params('freq') = freq;


%%
grid_h = 1:12;
% Sample period
for period = sample_tab.period(1)'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);

    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    % create another column for Rf with a standardized name
    idx_series_Rf = ismember({data.series}, series_Rf);
    T.Rf = T{:, data(idx_series_Rf).VariableNames}/100;
    if strcmp(freq, 'M')
        % pass, RF_FF is given in monthly terms
    elseif strcmp(freq, 'Q')
        T.Rf = T.Rf/4; % TBill given in annual terms
    elseif strcmp(freq, 'A')
        % pass, TBill given in annual terms
    else
        error('freq must be M, Q or A')
    end

    
    T.Rf_lag = lag(T.Rf, 1);
    T.Rf_lead = lead(T.Rf, 1);
    
    % compute inflation
    idx_series_CPI = ismember({data.series}, series_CPI);
    T.Pi = [NaN;...
        T{2:end,   data(idx_series_CPI).VariableNames}./...
        T{1:end-1, data(idx_series_CPI).VariableNames}];
    
    
    if log_approx
        T.Pi = log(T.Pi)*pcent;
    else
        T.Pi = (T.Pi -1)*pcent;
    end

    % reduce data to sample period
    T_ = T(idx_sample, :);

    summary_stats = [...
        [arrayfun(@(h) autocorr_(T_.Pi, h), grid_h), mean(T_.Pi), std(T_.Pi)];...
        [arrayfun(@(h) autocorr_(T_.Rf, h), grid_h), mean(T_.Rf), std(T_.Rf)]...
        ];
    array2table(summary_stats, 'VariableNames',...
        [strcat('rho_', arrayfun(@num2str, grid_h, 'UniformOutput', false)),...
        {'mean', 'std'}]);

    %{
    % check that inflation and Tbill are same frequency units
    [T_.Pi, T_.Rf]
    clf
    hold on
    plot(T_.date, T_.Rf*100)
    plot(T_.date, T_.Pi*100)
    return
    %}
    
    results_table = table;

    for t = 25:size(T_, 1) % minimum 24 observations for estimation

        Y = diff(T_.Pi(1:t-1)); % depdendent variable: first diff of inflation
        X = diff(T_.Rf(1:t-1)); % indepdendent variable: first diff of TBill
        
        N = size(X, 1);
        K = size(X, 2);

        % Functions for FGLS estimation
        % ****************************
        H_hat = @(theta) get_H_movavg(theta, N);
        Omega_inv_hat_ = @(theta) H_hat(theta)'*H_hat(theta);
    
        beta_FGLS_fun = @(theta) ((X'*Omega_inv_hat_(theta)*X)\eye(K))*...
            X'*Omega_inv_hat_(theta)*Y;
        sigma2_fun = @(beta) (Y-X*beta)'*(Y-X*beta)/N;
        Omega_inv_hat = @(theta) Omega_inv_hat_(theta)/sigma2_fun(beta_FGLS_fun(theta));
    
        %beta_FGLS = @(theta) X'*Omega_inv_hat(theta)*Y;
        % observe that errors are weighted by H_hat
    
        err = @(theta) (Y-X*beta_FGLS_fun(theta));
        errw = @(theta) H_hat(theta)*(Y-X*beta_FGLS_fun(theta));
        J = @(theta) sqrt(sum(err(theta).^2, 'omitnan')/N);
    
        % don't minimize SSR for theta, that is inconsistent. 
        % Implement maximum likelihood with FGLS, using *log*-likelihood
    
        L_fun = @(theta, beta)...
            -N/2*log(2*pi)-log(...
            sqrt(det(get_Omega_movavg(theta, sigma2_fun(beta), N))))-...
            (Y-X*beta)'*Omega_inv_hat(theta)*(Y-X*beta)/2/sigma2_fun(beta);
    
        L_fun_FGLS = @(theta)...
            -N/2*log(2*pi)-.5*sum(log(...
            eig(get_Omega_movavg(theta, sigma2_fun(beta_FGLS_fun(theta)), N))))-...
            (Y-X*beta_FGLS_fun(theta))'*...
            Omega_inv_hat(theta)*...
            (Y-X*beta_FGLS_fun(theta))/2;
    
        step = 1e-2; % precision parameter for theta_estimation
        
        % initial search
        theta_grid = -1.5:step:-.2;
        res_L = NaN(numel(theta_grid), 1);
        for theta = theta_grid
            clc
            fprintf('Computing log-likelihood function, %.2f%%',...
                find(theta_grid == theta)/numel(theta_grid)*100)
            
            res_L(theta_grid==theta) = L_fun_FGLS(theta);
        end
        
    
        subplot(1,2,2)
        plot(theta_grid, res_L)
        vline(theta_grid(res_L==max(res_L)))
        xlabel('\theta')
        ylabel('log-Likelihood function')
        drawnow
    
        step = 1e-3; % precision parameter for theta_estimation
        % initial search
        theta = theta_grid(res_L == max(res_L));
        theta_grid = theta-step*25:step:theta+25*step;
        res_L = NaN(numel(theta_grid), 1);
    
        for theta = theta_grid
            clc
            fprintf('Computing log-likelihood function, %.2f%%',...
                find(theta_grid == theta)/numel(theta_grid)*100)
    
            res_L(theta_grid==theta) = L_fun_FGLS(theta);
        end
    
        %clf
        plot(theta_grid, res_L)
        vline(theta_grid(res_L==max(res_L)))
        xlabel('\theta')
        ylabel('log-Likelihood function')
        drawnow
        
    
        step = 1e-4; % precision parameter for theta_estimation
        % initial search
        theta = theta_grid(res_L == max(res_L));
        theta_grid = theta-step*25:step:theta+25*step;
        res_L = NaN(numel(theta_grid), 1);
    
        for theta = theta_grid
            clc
            fprintf('Computing log-likelihood function, %.2f%%',...
                find(theta_grid == theta)/numel(theta_grid)*100)
    
            res_L(theta_grid==theta) = L_fun_FGLS(theta);
        end
    
        %clf
        plot(theta_grid, res_L)
        vline(theta_grid(res_L==max(res_L)))
        xlabel('\theta')
        ylabel('log-Likelihood function')
        drawnow
    
        beta_FGLS = beta_FGLS_fun(theta);
        
        % compare estimators: OLS and FGLS
        [beta_OLS, ~, ~, ~, STATS] = regress(Y, X);
        R_sq_OLS = STATS(1);
    
        Y_hat_FGLS = X*beta_FGLS_fun(theta);
        Y_hat_OLS = X*beta_OLS;
        Y_hat_1 = X * 1.00;

    
        err_OLS = Y-Y_hat_OLS;
        err_FGLS = Y-Y_hat_FGLS;
        err_1 = Y-Y_hat_1;
        theta_1 = autocorr_(err_1, 1);

        errw_FGLS = H_hat(theta)*(Y-Y_hat_FGLS); % rotated errors: FGLS
        errw_1 = H_hat(theta_1)*(Y-Y_hat_1); % rotated errors: beta=1
    
       
        % err_FGLS_t = u_t - theta * u_{t-1}
        u_0 = 0;
        u_FGLS = errw_FGLS + theta * [u_0; errw_FGLS(1:end-1)];
        fname = sprintf('residulas_FGLS_%s.xlsx', T_.date(t-1));
        fpath_FGLS = fullfile(root, 'FGLS', fname);
        writetable([T_(2:t-1, 'date'), array2table(err_FGLS), array2table(u_FGLS)],...
            fpath_FGLS)
        
    
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
        z1 = err_1; % set the observed signal as err_FGLS
        
        % variances
        % ***********
        % note: using unconditional var_z, var_v and var_eta causes
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
    
        var_z = var(z); % sample variance of signal, unconditional
        %var_z = arrayfun(@(t) var(z(t:t+periods_ma)), 1:N-periods_ma)'; % sample variance of signal, conditional
        %var_z = [var_z(1)*ones(periods_ma, 1); var_z];
    
        
        var_z1 = var(z1); % sample variance of signal, unconditional
        %var_z1 = arrayfun(@(t) var(z1(t:t+periods_ma)), 1:N-periods_ma)'; % sample variance of signal, conditional
        %var_z1 = [var_z1(1)*ones(periods_ma, 1); var_z1];
    
        autocorr_z = autocorr_(z, 1);
        %autocorr_z = arrayfun(@(t) autocorr_(z(t:t+periods_ma), 1), 1:N-periods_ma)';
        %autocorr_z = [autocorr_z(1)*ones(periods_ma, 1); autocorr_z];
    
        autocorr_z1 = autocorr_(z1, 1);
        %autocorr_z1 = arrayfun(@(t) autocorr_(z1(t:t+periods_ma), 1), 1:N-periods_ma)';
        %autocorr_z1 = [autocorr_z1(1)*ones(periods_ma, 1); autocorr_z1];
    
        var_eta = - autocorr_z .* var_z;
        var_eta1 = - autocorr_z1 .* var_z1;
    
        var_v = var_z - 2 * var_eta;
        var_v1 = var_z1 - 2 * var_eta1;
    
        if var_v < 0
            warning('Sample autocorrelation of signal lower than -0.50. Setting var_v to 10% of var_z')
            var_v = 0.1*var_z;
            var_v1 = 0.1*var_z1;
            var_eta = (var_z-var_v)/2;
            var_eta1 = (var_z1-var_v1)/2;
    
        elseif var_eta < 0
            warning('Sample autocorrelation of signal greater than 0.00. Setting var_v to 90% of var_z')
            var_v = 0.9*var_z;
            var_v1 = 0.9*var_z1;
            var_eta = (var_z-var_v)/2;
            var_eta1 = (var_z1-var_v1)/2;
    
        end
    
        
        % covariances
        % ***********
        Sigma_eta = eye(N) .* var_eta;
        Sigma_eta1 = eye(N) .* var_eta1;
        
        % construct Sigma_z
        Sigma_z_ = 3 * eye(N) - ones(N); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
        Sigma_z_ = Sigma_z_-spdiags(zeros(N, 3), -1:1, Sigma_z_);
        Sigma_z = Sigma_z_ .* var_eta + eye(N) .* var_v;
    
        % construct Sigma_z1
        Sigma_z1_ = 3 * eye(N) - ones(N); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
        Sigma_z1_ = Sigma_z1_-spdiags(zeros(N, 3), -1:1, Sigma_z1_);
        Sigma_z1 = Sigma_z1_ .* var_eta1 + eye(N) .* var_v1;
    
        % construct Sigma_eta_z
        Sigma_eta_z = 2 * eye(N) - ones(N); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
        Sigma_eta_z = Sigma_eta_z-spdiags(zeros(N, 2), 0:1, Sigma_eta_z);
        Sigma_eta_z = Sigma_eta_z .* var_eta;
    
        % construct Sigma_eta_z1
        Sigma_eta_z1 = 2 * eye(N) - ones(N); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
        Sigma_eta_z1 = Sigma_eta_z1-spdiags(zeros(N, 2), 0:1, Sigma_eta_z1);
        Sigma_eta_z1 = Sigma_eta_z1 .* var_eta1;
    
    
        if var_z ~= 1
            % if I use varying var_v and var_eta, can no long invert using this
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
            [~, inv_Sigma_z_] = get_H_movavg(-1, sample_size);
            % second component, simple inverse of identity matrix...
            
            % combining components following Miller, 1981 (Mathematics Magazine)
            inv_Sigma_z = inv_sum(...
                Sigma_z_*var_eta, eye(sample_size)*var_v,...
                inv_Sigma_z_/var_eta, eye(sample_size)/var_v...
                );
    
        end
        
        eta_0 = 0;
        E_eta_given_z = Sigma_eta_z * inv_Sigma_z * (z-0);
        Delta_eta = [E_eta_given_z(1)-eta_0; diff(E_eta_given_z)];
        Ev = z-Delta_eta;
    
        E_eta_given_z1 = Sigma_eta_z1 * inv_Sigma_z1 * (z1-0);
        Delta_eta1 = [E_eta_given_z1(1)-eta_0; diff(E_eta_given_z1)];
        Ev1 = z1-Delta_eta1;
    
    %{
    % sum of changes in alpha_t = -E real rate
    clf
    plot(T_.date(2:t), -cumsum(Ev)*100)
    return
    %}

    % total residuals of original model, correcting residulas using the sum of Ev
    resid = T_.Pi(2:t-1)-beta_FGLS*T_.Rf(2:t-1) - cumsum(Ev);
    
    alpha = cumsum(Ev) + mean(resid); % correct alpha by mean of resid, alpha_{t-1}
    resid_1 = T_.Pi(2:t-1)-1.0*T_.Rf(2:t-1) - cumsum(Ev); % replacing beta_FGLS with 1.0
    alpha_1 = cumsum(Ev) + mean(resid_1);

    EITB = alpha + beta_FGLS * T_.Rf(2:t-1); % in sample
    UITB = T_.Pi(2:t-1)-EITB; % in sample
    EITB_t = alpha + beta_FGLS * T_.Rf(t); % out of sample
    UITB_t = T_.Pi(t) - EITB_t; % out of sample

    EITB_1 = alpha_1 + 1.0 * T_.Rf(2:t-1);
    UITB_1 = T_.Pi(2:t-1)-EITB_1;
    EITB_t_1 = alpha_1 + 1.0 * T_.Rf(t);
    UITB_t_1 = T_.Pi(t) - EITB_t_1;
    
    SSR = sum(UITB.^2);
    SSR_1 = sum(UITB_1.^2);
    SS_tot = sum((T_.Pi(2:end)-mean(T_.Pi(2:end))).^2);
    R_sq = 1-SSR/SS_tot;
    R_sq_1 = 1-SSR_1/SS_tot;

    %{
    EITB = [T_(2:t, 'date'),...
        array2table(EITB*100, 'VariableNames', {'EITB'})];
    UITB = [T_(2:t, 'date'),...
        array2table(UITB*100, 'VariableNames', {'UITB'})];
    %}

        
    VariableNames = {...
        'Pi_t', 'EITB_{t-1}', 'EITB_{t-1}_1', 'UITB_t', 'UITB_t_1',...
        'N', 'K', 'beta_FGLS', 'beta_OLS', 'R_sq', 'R_sq_1', 'R_sq_OLS'};
    results_table_t = [...
        T_.Pi(t), EITB_t(end), EITB_t_1(end), UITB_t(end), UITB_t_1(end),...
        N, K, beta_FGLS, beta_OLS, R_sq, R_sq_1, R_sq_OLS];
    results_table_t = array2table(results_table_t, 'VariableNames', VariableNames);
    results_table_t.date = T_.date(t);
    results_table_t = results_table_t(:, [end, 1:end-1]);
    results_table = [results_table; results_table_t];

    % plot inflation and EITB as window expands
    subplot(1,2,1)
    cla
    hold on
    plot(T_.date, T_.Pi)
    plot(results_table.date, results_table.('EITB_{t-1}'))

    drawnow

    end

    
end

%%

fpath_FGLS_results_table = fullfile(root, 'FGLS_OOS.xlsx');
writetable(results_table, fpath_FGLS_results_table)


% EITB_t is expectation of Pi_t based on info at t-1
% UITB_t is inflation surprise at t based on info at t-1


%%
clc
clf
subplot(1,2,1)
bar(arrayfun(@(h) autocorr_(results_table.UITB_t, h), grid_h))
subplot(1,2,2)
plot(results_table.date, results_table.UITB_t)


%%
% it looks like the model misses inflation by one month.
% I checked and the timing is perfect...
clf
hold on
plot(T_.date, T_.Pi*100)
plot(results_table.date, results_table.('EITB_{t-1}')*100)
%plot(results_table.date(1:end-1), results_table.EITB_t(2:end)*100)



%% both OOS

fpath = fullfile(matlab_dir, 'NRC', 'arima_101_usa.csv')
arima_101 = readtable(fpath);
arima_101 = renamevars(arima_101, 'Pi_hat', 'Pi_hat_101');

fpath = fullfile(root, 'FGLS_OOS.xlsx')
TB_model = readtable(fpath);
%TB_model = renamevars(TB_model, 'Pi_hat', 'Pi_hat_TB');
%TB_model = renamevars(TB_model, 'err', 'epsilon_TB');

TB_model
arima_101

% get common dates
common_dates = intersect(arima_101.date, TB_model.date);
idx_arima = ismember(arima_101.date, common_dates);
idx_TB = ismember(TB_model.date, common_dates);
idx_T = ismember(T.date, common_dates);

clf
hold on
plot(T.date(idx_T), T.Pi(idx_T), 'r-')
plot(arima_101.date(idx_arima), arima_101.Pi_hat_101(idx_arima), 'k-')
plot(TB_model.date(idx_TB), TB_model.EITB__t_1_(idx_TB)*100, 'k--')

clf
hold on
%plot(T.date(idx_T), T.Pi(idx_T), 'r-')
plot(arima_101.date(idx_arima), arima_101.epsilon_101(idx_arima), 'k-')
plot(TB_model.date(idx_TB), TB_model.UITB_t(idx_TB)*100, 'r--')


%legend({'\pi', 'EIAR', 'EITB'})
legend({'EIAR', 'EITB'})
ylabel('Percent')

std(arima_101.Pi)
sqrt(mean((arima_101.epsilon_101(idx_arima)).^2))
sqrt(mean((TB_model.UITB_t(idx_TB)*100).^2))


corr(TB_model.EITB__t_1_(idx_TB), arima_101.Pi_hat_101(idx_arima))
corr(TB_model.UITB_t(idx_TB), arima_101.epsilon_101(idx_arima))


%%

clf
subplot(1,2,1)
bar(arrayfun(@(h) autocorr_(TB_model.UITB_t(idx_TB), h), grid_h))
subplot(1,2,2)
bar(arrayfun(@(h) autocorr_(arima_101.epsilon_101(idx_arima), h), grid_h))


%%

clf
clf
hold on
h=120
%plot(T.date(idx_T), T.Pi(idx_T), 'r-')
plot(arima_101.date(idx_arima), movmean(arima_101.epsilon_101(idx_arima), h), 'r-')
plot(TB_model.date(idx_TB), movmean(TB_model.UITB_t(idx_TB)*100, h), 'k--')















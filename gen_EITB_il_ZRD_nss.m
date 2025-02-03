
%{
In-sample procedure for TB model, FGLS estimation.
Data for Israel.

2024-06-11:
I've changed the estimation to a joint estimation of (theta, beta) via
Maximum likelihood (ML). For reference, I've kept computation of beta_FGLS
with the same method, but the ML is better.

The resulting beta is closer to -1 than FGLS. It also yields a low first
order autocorrelation for the unexpected inflation (UITB).

That it is still negative at the monthly frequency is beyond me.
This may mean that the changes in the short yields are driven by stronger
forces such as risk premia and/or liquidity factors.

it looks like all predictability is coming from change in the
expected real rate

2024-06-10:
I've added series_Rf, using the MAKAM yield curve. (MAKAM_yields)

Many times I get a newgative coefficient, which makes me think somehting is
wrong in this regression. theta does seem correct, and the resulting R_sq
etc shows this is something with the data itself.

Timing: I want to move to daily data, from the 16th to the 15th.
I need to make sure the timing i correct.

For example, MAKAM_yields15 is a monthly series, using the first business
day after the 15th of each month (say the 16th of Apr, 2024).
The inflation series dated 30th Apr regards Inflation between 1-30 Apr
and is announced on the 15th of May.

In order to match these I need to match these I need to push the MAKAM
yield to the end of April, which means not lagging inflation.

Lagging inflation adds one calendar month to inflation, so that Inflation
between 1-30 Apr, dated 30-Apr, will be dated 31-May.

Whatever I do, the prediction regression in Israel at monthly freq
produces a negative beta_FGLS. This has indeed has some prediction, but
the value of it is very low (I've made errors by including the t-residual
and effectively got exploded R_sq). The regression at annual freq
generates a coefficient of 1. The regression at quarterly freq
generates a coef of .25.

No matter if I use raw MAKAM, or the yield curve.

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
lib_data = fullfile('/media', username, 'D', 'data', 'BOI');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')


lib_israel = fullfile(lib_data, 'Israel');

lib_cpi = fullfile(lib_israel);
lib_ms = fullfile(lib_israel, 'monagg');
lib_mbase = fullfile(lib_israel, 'money_base');
lib_reserves = fullfile(lib_israel, 'reserves');

lib_gnp = fullfile(lib_israel, 'GNP');  % gross national product

lib_ea = fullfile(lib_israel, 'indprod');  % economic activity

lib_ir = fullfile(lib_israel);
lib_rr = fullfile(lib_israel, 'TSB_ZRD_nss', 'M');
%lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M');
lib_pop = fullfile(lib_israel, 'POP');
lib_consumption = fullfile(lib_israel, 'C');
%lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
lag_inflation = 1;
log_approx = false;
pcent = 100; % set as 100 if want percent
% keep pcent as 100, otherwise computation limits on log-likelihood
freq = 'M'; % frequency to be used

incl_xr = false;
incl_trgt = false;

series_CPI = 'CP_SA.M';
series_xr = 'RER_USD_ILS.D';
%series_ms = 'ZA108.M'; % waiting for Elad from MOS
series_ms = {'A346.M_E', 'ZA215.M'}; % ask Elad about RESNALNS, 
series_ea = 'CLS11.TPR_C_SA.Q_CHAINED'; % excl diamonds
series_trgt_min = 'INF_MIN_TRGT.D';
series_trgt_max = 'INF_MAX_TRGT.D';
series_rr = 'TSB_ZRD_001M'; % real rate
series_rr2 = 'TSB_ZRD_002M'; % real rate
series_rr3 = 'TSB_ZRD_003M'; % real rate


% total investment, incl financial sector
series_cx = 'AM_INV_TOTALT_Q_N';

%series_gnp = 'GNP.Q_N';
series_gnp = 'GNP.Q_FP';

series_ns = {'AM_NSTOCK_TOTALT_Q_N'};
series_C = {'C_@DUR.Q_FP_SA'}; % excl durable, should contain services

% not needed
%series_PD = {'PCEPINDG', 'PCEPIS'}; % I have the real series, don't really
%need this. (implicit deflator)
%series_pop = 'DEM_POP_AVE_CHN.Q';


incl_xr = false;
incl_trgt = false;


sample_tab = table(...
    [1:6]',...
    [datetime(1996,3,31), datetime(1996,3,31), datetime(2009,1,31),...
    datetime(1996,3,31), datetime(2000,1,31), datetime(2000,1,31),...
     ]',...
    [datetime(2019,12,31), datetime(2008,12,31), datetime(2019,12,31),...
    datetime(2024,4,30), datetime(2024,4,30), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});



%series_Rf = 'MAKAM_yields';
series_Rfs = {'MAKAM_raw', 'MAKAM_raw_15', 'MAKAM_raw_15mean',...
        'MAKAM_yields', 'MAKAM_yields_15', 'MAKAM_yields_15mean',...
        'TELBOR', 'TELBOR_15', 'TELBOR_15mean',...
        };
for series_Rf = series_Rfs

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});
if strcmp(series_Rf, 'MAKAM_raw')
    % X-month MAKAM, annualized average of calendar month.
    % raw data gets the closest bond to X-months, which inserts noise, due
    % to liquidity as well as tick size. Not recommended.

    series_Rf = 'TSB_BAGR_MAKAM_01M.M'; % 1 month MAKAM, annualized
    %series_Rf = 'TSB_BAGR_MAKAM_02M.M'; % 2 month MAKAM, annualized
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized

    lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M', 'with_metadata');
    
elseif strcmp(series_Rf, 'MAKAM_raw_15')
    % X-month MAKAM, annualized from daily (not average) after 15 of the month
    % raw data gets the closest bond to X-months, which inserts noise, due
    % to liquidity as well as tick size. Not recommended.

    series_Rf = 'TSB_BAGR_MAKAM_01M.M_15'; % 1 month
    %series_Rf = 'TSB_BAGR_MAKAM_02M.M_15'; % 2 month
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M_15'; % 3 month

    lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM_15', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_raw_15mean')
    % X-month MAKAM, annualized from daily (average) after 15 of the month
    % to the 15th of the next month.
    % raw data gets the closest bond to X-months, which inserts noise, due
    % to liquidity as well as tick size. Not recommended.

    series_Rf = 'TSB_BAGR_MAKAM_01M.M_15mean'; % 1 month
    %series_Rf = 'TSB_BAGR_MAKAM_02M.M_15mean'; % 2 month
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M_15mean'; % 3 month

    lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM_15mean', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_yields')
    % X-month MAKAM, annualized average of calendar month --
    % yield curve estimation using NS model.

    series_Rf = 'MAKAM_yields_M01.M'; % 1 month
    lib_ir = fullfile(lib_israel, 'MAKAM_yields', 'M', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_yields_15')
    % X-month MAKAM, annualized from daily (not average) after 15 of the
    % month --
    % yield curve estimation using NS model.

    series_Rf = 'MAKAM_yields_M01.M_15'; % 1 month
    lib_ir = fullfile(lib_israel, 'MAKAM_yields_15', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_yields_15mean')
    % X-month MAKAM, annualized from daily (average) after 15 of the month
    % to the 15th of the next month --
    % yield curve estimation using NS model.

    series_Rf = 'MAKAM_yields_M01.M_15mean'; % 1 month
    lib_ir = fullfile(lib_israel, 'MAKAM_yields_15mean', 'with_metadata');

elseif strcmp(series_Rf, 'TELBOR')
    % X-month MAKAM, annualized average of calendar month --
    % TEBOR rates

    series_Rf = 'BL.TELBOR_01M.M'; % 1 month
    lib_ir = fullfile(lib_israel, 'TELBOR', 'M', 'with_metadata');

elseif strcmp(series_Rf, 'TELBOR_15')
    % X-month MAKAM, annualized from daily (not average) after 15 of the
    % month --
    % yield curve estimation using NS model.

    series_Rf = 'BL.TELBOR_01M.M_15'; % 1 month
    lib_ir = fullfile(lib_israel, 'TELBOR_15', 'with_metadata');
    
elseif strcmp(series_Rf, 'TELBOR_15mean')
    % X-month MAKAM, annualized from daily (average) after 15 of the month
    % to the 15th of the next month --
    % yield curve estimation using NS model.

    series_Rf = 'BL.TELBOR_01M.M_15mean'; % 1 month
    lib_ir = fullfile(lib_israel, 'TELBOR_15mean', 'with_metadata');

else
    error('Unknown %s', series_Rf)
end


series_Rf2 = strrep(series_Rf, '01', '02');
series_Rf3 = strrep(series_Rf, '01', '03');

if strcmp(freq, 'M')
    % pass
elseif strcmp(freq, 'Q')
    series_Rf = strrep(series_Rf, '01', '03');
    series_rr = strrep(series_rr, '01', '03');
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M_15';
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized
elseif strcmp(freq, 'A')
    series_Rf = strrep(series_Rf, '01', '12');
    series_rr = strrep(series_rr, '01', '12');
    %series_Rf = 'TSB_BAGR_MAKAM_12M.M'; % 12 month MAKAM, annualized
else
    error('freq must be M, Q or A')
end


% src: BI, Bank of Israel old website; FAME
data = struct(...
    'series', {...
    'CP_NSA.M',... CPI, not seasonally adjusted, BLS
    'CP_SA.M',... Price index of all urban consumers
    series_rr,...
    series_rr2,...
    series_rr3,...
    series_Rf,...
    series_Rf2,...
    series_Rf3,...
    series_xr,...
    'ZA108.M',...
    'A346.M_E',...
    'ZA215.M',...
    series_ea,...
    series_trgt_min,...
    series_trgt_max,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_rr,...
    lib_rr,...
    lib_rr,...
    lib_ir,...
    lib_ir,...
    lib_ir,...
    lib_israel,...
    lib_mbase,...
    lib_mbase,...
    lib_reserves,...
    lib_ea,...
    lib_israel,...
    lib_israel,...
    %'',...
    },...
    'src', {...
    'FAME',...
    'FAME',...
    'NONE',...
    'NONE',...
    'NONE',...
    'FAME',...
    'FAME',...
    'FAME',...
    'FAME',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'FAME',...
    'FAME',...
    %'',...
    },...
    'VariableNames', cell(1,15)...
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
    [T_, metadata_] = load_data_files(data(i).series,...
        lib_data_series, data(i).src);
    data(i).VariableNames = T_.Properties.VariableNames(2);

    if strcmp(data(i).series, series_Rf) ||...
            strcmp(data(i).series, series_Rf2) ||...
            strcmp(data(i).series, series_Rf3)
        T_.date.Day = 1;
        T_.date = T_.date+calmonths(1);
    end
    if strcmp(data(i).series, series_rr) ||...
            strcmp(data(i).series, series_rr2) ||...
            strcmp(data(i).series, series_rr3)
        T_.date.Day = 1;
        T_.date = T_.date+calmonths(1);
    end
    if strcmp(data(i).series, series_xr)
        T_ = table2timetable(T_);
        T_ = retime(T_, 'daily', 'spline');
        T_ = timetable2table(T_);
        %T_{:, end} = fillmissing(T_{:, end}, 'previous');
    end

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


%
% error in BOI CP_SA series... in Jan-2000.
% seasonally adjust CP_NSA for period not available from CBS, adding 12
% months for sa-estimation
idx_CP_SA = 1:find(~isnan(T.CP_SA0x2EM), 1, 'first')+12;
CP_SA = sa_adj(T.CP_NSA0x2EM(idx_CP_SA), 12);
% set CP_SA as the official CBS series of the CPI
T.CP_SA = T.CP_SA0x2EM; % error in the official BOI series
%}


%idx_CP_SA = T.date.Year > 1990;
%CP_SA = sa_adj(T.CP_NSA0x2EM, 12);
%T.CP_SA(idx_CP_SA) = CP_SA;
%}
%
% complete this series with the above CP_SA
T.CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1) = ...
    CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1);
%}



grid_h = 1:12;
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
    grid_h = 1:4;
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

%T = fillmissing(T, 'previous');



params = containers.Map();
params('constant') = true;
params('weight') = false;
params('freq') = freq;


% Sample period
for period = sample_tab.period(1:end)'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);
    vnames_ML_FGLS = {'theta', 'beta_rf'};

    fname_ML = 'ML_';
    fname_FGLS = 'FGLS_';

    lib_results = fullfile(root, 'ML_FGLS_il_ZRD_nss', 'Rf');

    if incl_xr
        lib_results = [lib_results, '_xr'];
        fname_ML = [fname_ML, 'xr_'];
        fname_FGLS = [fname_FGLS, 'xr_'];
        vnames_ML_FGLS = [vnames_ML_FGLS, {'beta_xr'}];
    end
    if incl_trgt
        lib_results = [lib_results, '_trgt'];        
        fname_ML = [fname_ML, 'trgt_'];
        fname_FGLS = [fname_FGLS, 'trgt_'];
        vnames_ML_FGLS = [vnames_ML_FGLS, {'beta_trgt'}];
    end

    vnames_ML_FGLS = [vnames_ML_FGLS, {'mean_Pi', 'mean_Rf', 'rho1_Pi', 'rho1_Rf'}];

    fname_ML = [fname_ML, series_Rf];
    fname_FGLS = [fname_FGLS, series_Rf];

    fname_ML = sprintf('%s_%s--%s.xlsx', fname_ML, min_date, max_date);
    fname_FGLS = sprintf('%s_%s--%s.xlsx', fname_FGLS, min_date, max_date);
    
    fpath_ML = fullfile(lib_results, fname_ML);
    fpath_FGLS = fullfile(lib_results, fname_FGLS);

    % data to save for R to run a KalmanFilter procedure
    fname_R_in = sprintf('%s--%s_residulas_ML_FGLS.xlsx', min_date, max_date);
    fpath_R_in = fullfile('/home' username, 'R', 'Fama_Gibbons_1982', 'data', 'il', 'ZRD_nss', series_Rf, fname_R_in);

    if exist(fpath_R_in, 'file') == 2
        continue
    end

    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    idx = ismember({data.series}, {series_trgt_min, series_trgt_max});
    T.trgt = mean(T{:, [data(idx).VariableNames]}, 2, 'omitnan');

    idx = ismember({data.series}, series_Rf);
    T.Rf = T{:, data(idx).VariableNames}/100;
    idx = ismember({data.series}, series_rr);
    T.rr = T{:, data(idx).VariableNames}/100;
    

    if strcmp(freq, 'M')
        T.Rf = T.Rf/12; % data in annual terms -- to monthly
        T.rr = T.rr/12; % data in annual terms -- to monthly
        T.trgt = T.trgt/12;
    elseif strcmp(freq, 'Q')
        T.Rf = T.Rf/4;
        T.trgt = T.trgt/4;
        T.rr = T.rr/4;
    elseif strcmp(freq, 'A')
        % pass
    else
        error('freq must be M, Q or A')
    end
        
    T.Pi = [NaN;...
        T.CP_SA(2:end, :)./...
        T.CP_SA(1:end-1, :)];

    %{
    idx = ismember({data.series}, series_CPI);
    T.Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    %}
    idx = ismember({data.series}, series_xr);
    T.xr = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];

    
    if log_approx
        T.Pi = log(T.Pi)*pcent;
        T.xr = log(T.xr)*pcent;
    else
        T.Pi = (T.Pi -1)*pcent;
        T.xr = (T.xr -1)*pcent;
    end
    T.Rf = T.Rf*pcent;
    T.rr = T.rr*pcent;


    T_ = T(idx_sample, :);

    summary_stats = [...
        [arrayfun(@(h) autocorr_(T_.Pi, h), grid_h), mean(T_.Pi), std(T_.Pi)];...
        [arrayfun(@(h) autocorr_(T_.Rf, h), grid_h), mean(T_.Rf), std(T_.Rf)]...
        ];
    array2table(summary_stats, 'VariableNames',...
        [strcat('rho_', arrayfun(@num2str, grid_h, 'UniformOutput', false)),...
        {'mean', 'std'}])

    
    % regular model
    Y = T_.Pi;
    X = T_{:, {'Rf'}}-T_.rr;
    % model in first difference
    Y = diff(Y, 1);
    X = diff(X, 1);

    if incl_xr
        X = [X, T_{:, {'xr'}}];
    end
    if incl_trgt
        X = [X, T_{:, {'trgt'}}];
    end

    if any(isnan([Y(:); X(:)]))
        continue
        error('Missing values.')
    end
    
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
    [min_date, max_date]
    series_Rf
    % init conditions
    theta_prev = -1;
    theta_ML = -.8;
    beta_ML = 1;

    precision = 1e-5;
    precision_theta = 1e-2;
    precision_beta = 1e-2;
    spectrum_theta = [-1, -.2];
    spectrum_beta = [-2, 2];
    
    counter = 0;
    while abs(theta_ML-theta_prev) > precision
        counter = counter + 1;
        theta_prev = theta_ML;
        
        [beta_ML, ~] = my_fminunc(@(beta_) -L_fun(theta_ML, beta_), beta_ML, spectrum_beta, precision_beta);
        [theta_ML, ~] = my_fminunc(@(theta_) -L_fun(theta_, beta_ML), theta_ML, spectrum_theta, precision_theta);
        
        if mod(counter, 5) == 0
            precision_theta = max(precision_theta/10, 1e-5);
            precision_beta = max(precision_beta/10, 1e-5);
        end

        spectrum_beta = [beta_ML-10*precision_beta, beta_ML+10*precision_beta];
        spectrum_theta = [theta_ML-10*precision_theta, theta_ML+10*precision_theta];
        spectrum_theta = max(spectrum_theta, -1);

        fprintf('Fixed point estimation: theta_{ML}=%.4f  beta_{ML}=%.4f\n',...
            theta_ML, beta_ML)
    end
    
    % Maximum likeloihhod estimation, on theta + FGLS assumption; iteratively
    % *******************************************************************
    
    % init conditions
    theta_prev = -1;
    theta_FGLS = -.8;
    precision = 1e-5;
    precision_theta = 1e-2;
    spectrum_theta = [-2.5, -.2];
    
    counter = 0;
    while abs(theta_FGLS-theta_prev) > precision
        counter = counter + 1;
        theta_prev = theta_FGLS;
        
        [theta_FGLS, ~] = my_fminunc(@(theta_) -L_fun_FGLS(theta_), theta_FGLS, spectrum_theta, precision_theta);
        beta_FGLS = beta_FGLS_fun(theta_FGLS);
        
        if mod(counter, 5) == 0
            precision_theta = max(precision_theta/10, 1e-5);
        end
        spectrum_theta = [theta_FGLS-10*precision_theta, theta_FGLS+10*precision_theta];
        spectrum_theta = max(spectrum_theta, -1);
        
        fprintf('Fixed point estimation: theta_{FGLS}=%.4f  beta_{FGLS}=%.4f\n',...
            theta_FGLS, beta_FGLS)

    end

    
    %{
    %% plot the joint log-likelihood function

    %min_theta = min(theta_ML, theta_FGLS)-abs(theta_ML-theta_FGLS)*5;
    %max_theta = max(theta_ML, theta_FGLS)+abs(theta_ML-theta_FGLS)*5;
    %min_beta = min(beta_ML, beta_FGLS)-abs(beta_ML-beta_FGLS)*5;
    %max_beta = max(beta_ML, beta_FGLS)+abs(beta_ML-beta_FGLS)*5;

    min_theta = min(-1, min(theta_FGLS, theta_ML)*1.1);
    max_theta = -.2;
    min_beta = min(-1.5, min(beta_FGLS(1), beta_ML(1))-.1);
    max_beta = 1.5;

    grid_theta = linspace(min_theta, max_theta, 100);
    grid_beta = linspace(min_beta, max_beta, 100);

    grid_theta2 = repmat(grid_theta, numel(grid_beta), 1)'; % ax=1
    grid_beta2 = repmat(grid_beta, numel(grid_theta), 1); % ax=2
    l = arrayfun(@(i) L_fun(grid_theta2(i), grid_beta2(i)), 1:numel(grid_theta2));
    l = reshape(l, [numel(grid_theta), numel(grid_beta)]);
    
    % winsorize l
    l(l<prctile(l(:), 30)) = NaN;
    
    % check reshape
    %all(all(reshape(grid_theta2(:), [numel(grid_theta), numel(grid_beta)])==grid_theta2))
    %all(all(reshape(grid_beta2(:), [numel(grid_theta), numel(grid_beta)])==grid_beta2))
    
    clf
    hold on
    contourf(grid_theta, grid_beta, l')
    colorbar
    scatter(theta_ML, beta_ML(1), 'r*')
    scatter(theta_FGLS, beta_FGLS(1), 'k*')
    xlabel('\theta (moving average parameter)')
    ylabel('\beta')
    title('Joint log-likelihood function',...
        sprintf('Israel, %s -- %s', min_date, max_date))
    legend({'', 'ML', 'FGLS'})
    fname_fig = strrep(fname_ML, 'ML_', '');
    fname_fig = strrep(fname_fig, '.xlsx', '.jpg');
    fpath_fig = fullfile(root, 'gfx', 'IL', 'log-likelihood_ZRD_nss', fname_fig);
    export_fig(fpath_fig)

    %continue
    %%
    %}

    %%

    %
    [beta_OLS, beta_OLS_INT, ~, ~, STATS] = regress(Y, X);
    %[beta_OLS_INT(:, 1), beta_OLS, beta_OLS_INT(:, 2)]

    % STATS(1) % R_sq
    Y_hat_ML = X*beta_ML;
    Y_hat_FGLS = X*beta_FGLS;
    Y_hat_OLS = X*beta_OLS;
    Y_hat_1 = X * 1.00;

    err_ML = Y-Y_hat_ML;
    err_FGLS = Y-Y_hat_FGLS;
    err_OLS = Y-Y_hat_OLS;
    err_1 = Y-Y_hat_1;

    errw_ML = H_hat(theta_ML)*(Y-Y_hat_ML);
    errw_FGLS = H_hat(theta_FGLS)*(Y-Y_hat_FGLS);

    SSR_ML_diff = sum(err_ML.^2);
    SSR_FGLS_diff = sum(err_FGLS.^2);
    SSR_1_diff = sum(err_1.^2);
    SS_tot_diff = sum((Y-mean(Y)).^2);

    R_sq_ML_diff = 1-SSR_ML_diff/SS_tot_diff;
    R_sq_FGLS_diff = 1-SSR_FGLS_diff/SS_tot_diff;
    R_sq_OLS_diff = STATS(1);
    R_sq_1_diff = 1-SSR_1_diff/SS_tot_diff;
    %{
    %% check  R_sq of diff-regression
    clc
    [min_date, max_date]
    [beta_ML, beta_FGLS, beta_OLS, 1]
    [R_sq_ML_diff, R_sq_FGLS_diff, R_sq_OLS_diff, R_sq_1_diff,]

    clf
    hold on
    plot(T_.date(2:end), err_FGLS)
    plot(T_.date(2:end), err_1)
    %%
    %}


    %{
    % check autocorrelation in residuals
    acf_err_ML = arrayfun(@(h) autocorr_(err_ML, h), grid_h);
    acf_err_FGLS = arrayfun(@(h) autocorr_(err_FGLS, h), grid_h);
    acf_err_OLS = arrayfun(@(h) autocorr_(err_OLS, h), grid_h);
    acf_errw_ML = arrayfun(@(h) autocorr_(errw_ML, h), grid_h);
    acf_errw_FGLS = arrayfun(@(h) autocorr_(errw_FGLS, h), grid_h);
    acf_err_1 = arrayfun(@(h) autocorr_(err_1, h), grid_h);
    clf
    %bar(grid_h, [acf_u_FGLS; acf_u_1; acf_err_OLS; acf_err_FGLS])
    %legend({'u_{FGLS}', 'u_{FAMA}', 'err_{OLS}', 'err_{FGLS}'})
    bar(grid_h, [acf_err_FGLS; acf_errw_FGLS])
    legend({'err_{FGLS}', 'errw_{FGLS}'})
    %bar(grid_h, [acf_err_1; acf_errw_1])
    %legend({'err_{1}', 'errw_{1}'})
    %bar(grid_h, acf_err_OLS)
    %legend({'err_{OLS}'})
    % here we see that Fama's residulas, and the breakdown to u are very
    % close. The autocorrelation in u is still substantial.
    return
    %}
    
    
    if exist(fullfile('/home' username, 'R', 'Fama_Gibbons_1982', 'data', 'il', 'ZRD_nss', series_Rf), 'dir') ~= 2
        mkdir(fullfile('/home' username, 'R', 'Fama_Gibbons_1982', 'data', 'il', 'ZRD_nss', series_Rf))
    end

    writetable([T_(2:end, 'date'),...
        array2table(err_ML),...
        array2table(err_FGLS),...
        array2table(err_1)],...
        fpath_R_in)
    continue
    %%
    
    % output of KalmanFilter procedure from R
    fname_R_out = sprintf('%s--%s_KalmanFilter.csv', min_date, max_date);
    fpath_R_out = fullfile('/home', username, 'R', 'Fama_Gibbons_1982', 'results', 'il', 'ZRD_nss',...
        series_Rf, fname_R_out);

    if exist(fpath_R_out, 'file') ~= 2
        %pass
    else
        R_out = readtable(fpath_R_out);
        R_out.Properties.VariableNames = {'date', 'u_t', 'u_t_', 'v_t'};
        R_out.date = T_.date;
    end


    %{
    ***************************************************
    Wandering intercept procedure
    ***************************************************
    code extract from wandering_intercept.m, uses residulas_FGLS.xlsx
    

    
    z_t = u_t - u_{t-1} + v_t  # observed signal (residual)
    
    u_t ~ N(0, var_u)
    v_t ~ N(0, var_v)
    
    cov(u_t, u_s) = 0 for all s ~= t
    cov(v_t, v_s) = 0 for all s ~= t
    cov(u_t, v_s) = 0 for any s,t
    
    get conditional expectation of u given observed z
    % extract v (conditional expectation...)
    
    %}
        
    z = err_ML; % set the observed signal as ML
    z1 = err_1; % set the observed signal as 1
    
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

    var_u = - autocorr_z .* var_z;
    var_u1 = - autocorr_z1 .* var_z1;

    var_v = var_z - 2 * var_u;
    var_v1 = var_z1 - 2 * var_u1;

    if var_v < 0
        warning('autocorrelation of signal lower than -0.50. Setting var_v to 10% of var_z')
        var_v = 0.1*var_z;
        var_v1 = 0.1*var_z1;
        var_u = (var_z-var_v)/2;
        var_u1 = (var_z1-var_v1)/2;
    end

    %var_v = max(var_v, 0);
    %var_v1 = max(var_v1, 0);
    %var_u = (var_z-var_v)/2;
    %var_u1 = (var_z1-var_v1)/2;

    %autocorr_(z, 1) % sample 1st autocovariance
    %autocorr_(z, 1) * var_z % sample estimate of cov(z_t, z_t-1)
    %var_z/var_v % signal to noise ratio
    
    
    % covariances
    % ***********
    Sigma_u = eye(N) .* var_u;
    Sigma_u1 = eye(N) .* var_u1;
    
    % construct Sigma_z
    Sigma_z_ = 3 * eye(N) - ones(N); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
    Sigma_z_ = Sigma_z_-spdiags(zeros(N, 3), -1:1, Sigma_z_);
    Sigma_z = Sigma_z_ .* var_u + eye(N) .* var_v;

    % construct Sigma_z1
    Sigma_z1_ = 3 * eye(N) - ones(N); % [(2, -1, 0,..., 0); (0, 2, -1,...,0);...]
    Sigma_z1_ = Sigma_z1_-spdiags(zeros(N, 3), -1:1, Sigma_z1_);
    Sigma_z1 = Sigma_z1_ .* var_u1 + eye(N) .* var_v1;

    % construct Sigma_u_z
    Sigma_u_z = 2 * eye(N) - ones(N); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
    Sigma_u_z = Sigma_u_z-spdiags(zeros(N, 2), 0:1, Sigma_u_z);
    Sigma_u_z = Sigma_u_z .* var_u;

    % construct Sigma_u_z1
    Sigma_u_z1 = 2 * eye(N) - ones(N); % [(1, -1, 0,..., 0); (0, 1, -1,...,0);...]
    Sigma_u_z1 = Sigma_u_z1-spdiags(zeros(N, 2), 0:1, Sigma_u_z1);
    Sigma_u_z1 = Sigma_u_z1 .* var_u1;


    if numel(var_z) ~= 1
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
            Sigma_z_*var_u, eye(N)*var_v,...
            inv_Sigma_z_/var_u, eye(N)/var_v...
            );
        inv_Sigma_z1 = inv_sum(...
            Sigma_z1_*var_u1, eye(N)*var_v1,...
            inv_Sigma_z_/var_u1, eye(N)/var_v1...
            );

    end

    u_0 = 0;
    E_u_given_z = Sigma_u_z * inv_Sigma_z * (z-0);
    Delta_u = [E_u_given_z(1)-u_0; diff(E_u_given_z)];
    Ev = z-Delta_u;

    E_u_given_z1 = Sigma_u_z1 * inv_Sigma_z1 * (z1-0);
    Delta_u1 = [E_u_given_z1(1)-u_0; diff(E_u_given_z1)];
    Ev1 = z1-Delta_u1;
    %{
    % Comparison of KalmanFilter method with CondNorm method
    % ******************************************************
    %%
    if exist(fpath_R_out, 'file') == 2
        clc
        acf_err_OLS = arrayfun(@(h) autocorr_(err_OLS, h), grid_h);
        acf_err_FGLS = arrayfun(@(h) autocorr_(err_FGLS, h), grid_h);
        acf_errw_FGLS = arrayfun(@(h) autocorr_(errw_FGLS, h), grid_h);
        acf_E_u_given_z = arrayfun(@(h) autocorr_(E_u_given_z, h), grid_h);
        acf_R_out_u_t = arrayfun(@(h) autocorr_(R_out.u_t(2:end), h), grid_h);
    
        
        clf
        subplot(1,2,1)
        hold on
        plot(R_out.date(2:end), R_out.u_t(2:end), 'LineWidth', 2)
        plot(T_.date(2:end), E_u_given_z, '--')
        %title(sprintf('acf1: KalmanFilter=%.3f; CondNormal=%.3f',...
        %    autocorr_(R_out.u_t(2:end), 1), autocorr_(E_u_given_z, 1)))
        title(sprintf('mean: KalmanFilter=%.3f; CondNormal=%.3f',...
            mean(R_out.u_t(2:end)), mean(E_u_given_z)))
        legend({'u_t (KalmanFilter)', 'u_t (CondNormal)'})
    
        subplot(1,2,2)
        %bar(grid_h, [acf_u_FGLS; acf_u_FAMA; acf_err_OLS; acf_err_FGLS])
        %legend({'u_{FGLS}', 'u_{FAMA}', 'err_{OLS}', 'err_{FGLS}'})
        bar(grid_h, [acf_E_u_given_z; acf_R_out_u_t])
        legend({'u_t (KalmanFilter)', 'u_t (CondNormal)'})
        ylim([-1 1])
        fname_comparing_KalmanFilter2CondNormal = sprintf('IL__%s--%s.jpg',...
            min_date, max_date);
        fpath_comparing_KalmanFilter2CondNormal = fullfile('.', 'gfx',...
            'comparing_KalmanFilter2CondNormal', 'il', series_Rf,...
            fname_comparing_KalmanFilter2CondNormal);
        export_fig(fpath_comparing_KalmanFilter2CondNormal)
    end
    continue

    %%
    %}

    %
    clc; clf
    %plot(T_.date(3:end), -cumsum(Ev)*100)

    %plot(T_.date, T_.Pi*100)
    %resid = T_.Pi(2:end)-beta_FGLS*T_.Rf(2:end) - cumsum(Ev);
    %alpha = cumsum(Ev) + mean(resid);
    resid = T_.Pi-beta_ML*(T_{:, {'Rf'}}-T_.rr) - [0; cumsum(Ev)];
    alpha_t = [0; cumsum(Ev)] + mean(resid);

    %resid1 = T_.Pi(2:end)-1*T_.Rf(2:end) - cumsum(Ev1);
    %alpha1 = cumsum(Ev1) + mean(resid1);
    resid1 = T_.Pi-1*(T_{:, {'Rf'}}-T_.rr) - [0; cumsum(Ev1)];
    alpha1_t = [0; cumsum(Ev1)] + mean(resid1);

    EITB = alpha_t(1:end-1) + beta_FGLS * (T_{2:end, {'Rf'}}-T_.rr(2:end));
    UITB = T_.Pi(2:end)-EITB;
    
    EITB1 = alpha1_t(1:end-1) + 1.0 * (T_{2:end, {'Rf'}}-T_.rr(2:end));
    UITB1 = T_.Pi(2:end)-EITB1;
    %
    clf
    subplot(1,2,1)
    hold on
    plot(T_.date(2:end), EITB*100)
    plot(T_.date(2:end), EITB1*100)
    plot(T_.date, T_.Pi*100, 'k--')
    legend({'EITB', 'EITB1', '\pi', })
    
    subplot(1,2,2)
    hold on

    %plot(T_.date(3:end), UITB)
    %plot(T_.date(3:end), UITB1)
    h=0;
    plot(T_.date(3+h:end), arrayfun(@(t) mean(UITB(t:t+h)), 1:N-h-1)')
    h=12;
    plot(T_.date(3+h:end), arrayfun(@(t) mean(UITB(t:t+h)), 1:N-h-1)')
    title(sprintf('E(UITB)=%.3f', mean(UITB)))
    %}


    drawnow
    return

    SSR = sum(UITB.^2);
    SSR_1 = sum(UITB1.^2);
    SS_tot = sum((T_.Pi(3:end)-mean(T_.Pi(3:end))).^2);
    R_sq = 1-SSR/SS_tot;
    R_sq_OLS = STATS(1);
    R_sq_1 = 1-SSR_1/SS_tot;
    %}
    

    clf
    hold on
    %plot(T_.date(2:end), err_1)
    %plot(T_.date(2:end), err_FGLS)
    
    %plot(T_.date(3:end), EITB1)
    %plot(T_.date(3:end), EITB)
    %plot(T_.date, T_.Pi, 'k--')

    plot(T_.date(2:end), Y_hat_FGLS)
    plot(T_.date(2:end), Y_hat_1)
    plot(T_.date(2:end), Y, 'k--')

    %%
    
    EITB = [T_(2:end, 'date'),...
        array2table(EITB*100/pcent, 'VariableNames', {'EITB'})];
    UITB = [T_(2:end, 'date'),...
        array2table(UITB*100/pcent, 'VariableNames', {'UITB'})];
    %{
    EITB1 = [T_(2:end, 'date'),...
        array2table(EITB1*100/pcent, 'VariableNames', {'EITB1'})];
    UITB1 = [T_(2:end, 'date'),...
        array2table(UITB1*100/pcent, 'VariableNames', {'UITB1'})];
    %}
    
    lib_save_EITB = fullfile(root, 'EITB_il', series_Rf);
    if exist(lib_save_EITB, "dir") ~= 2
        mkdir(lib_save_EITB)
    end

    lib_save_UITB = fullfile(root, 'UITB_il', series_Rf);
    if exist(lib_save_UITB, "dir") ~= 2
        mkdir(lib_save_UITB)
    end

    fname = sprintf('EITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB = fullfile(lib_save_EITB, fname);
    writetable(EITB, fpath_EITB)
    
    fname = sprintf('UITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB = fullfile(lib_save_UITB, fname);
    writetable(UITB, fpath_UITB)

    %{
    fname = sprintf('EITB1__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB1 = fullfile(root, 'EITB_il', fname);
    writetable(EITB1, fpath_EITB1)
    
    fname = sprintf('UITB1__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB1 = fullfile(root, 'UITB_il', fname);
    writetable(UITB1, fpath_UITB1)
    %}
    
end
clear metadata T_ T
end


return

%%

period = 1;
min_date = sample_tab.min_date(period);
max_date = sample_tab.max_date(period);


fpath = fullfile(matlab_dir, 'NRC', 'arima_101_il.csv');
arima_101 = readtable(fpath);
arima_101 = renamevars(arima_101, 'Pi_hat', 'Pi_hat_101');
arima_101.date = arima_101.date+1;
arima_101.date = arima_101.date - calmonths(1);
arima_101.date = arima_101.date-1;

fname = sprintf('EITB__%s-%s.%s.csv', min_date, max_date, freq);
fpath_EITB = fullfile(root, 'EITB_il', fname);
fname = sprintf('UITB__%s-%s.%s.csv', min_date, max_date, freq);
fpath_UITB = fullfile(root, 'UITB_il', fname);

EITB = readtable(fpath_EITB);
UITB = readtable(fpath_UITB);

TB_model = join(EITB, UITB);

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





%{

Conditional model to predict next month inflation using one-month interest
rates. 


The OOS-RMSE is lower when limiting the sample to a few years because the
parameters allowed to change over time. I varied the sample size between
2-20 years, and RMSE is lowest around 5 years or less. Observe that the
OOS-R^2 is not necessarily good, which implies in certain periods interest
rates do not contain additional information above a constant. This is fine,
because in some years interest rates do not move by much.

pp. 549
(6) Pi_t = -b_0 - b_1 dA_t - b_2 dRf_t + b_3 dM_t + eta_t

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
%{
% extended dates
sample_tab = table(...
    [1:8]',...
    [datetime(1953,1,31), datetime(1953,1,31), datetime(2000,1,31),...
     datetime(1953,1,31), datetime(2007,7,31), datetime(1953,12,31),...
     datetime(1954,1,31), datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(2019,12,31), datetime(1977,12,31),...
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


freq = 'M'; % frequency to be used
%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

TB_model_results = [];


series_CPI = 'CPIAUCSL';
series_Rf = 'RF_FF';


% src: BI, Bank of Israel old website; FAME
data = struct(...
    'series', {...
    series_CPI,... Price index of all urban consumers
    series_Rf,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_ir,...
    %'',...
    },...
    'src', {...
    'FRED',...
    'FRED',...
    %'',...
    },...
    'VariableNames', cell(1,2)...
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



sample_tab = T.date(T.date.Year>=1953);
sample_tab = [repmat(datetime(1951,1,31), size(sample_tab)), sample_tab];
%sample_tab = [datetime(sample_tab.Year-5, sample_tab.Month, sample_tab.Day), sample_tab];
sample_tab = array2table(sample_tab, 'VariableNames', {'min_date', 'max_date'});
sample_tab.period = [1:size(sample_tab, 1)]';
sample_tab = sample_tab(:, {'period', 'min_date', 'max_date'});



% Sample period
for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);


    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames};
    Rf_lag = [NaN; T{1:end-1, data(idx).VariableNames}];
    Rf_lead = [T{2:end, data(idx).VariableNames}; NaN];

    
    
    idx = ismember({data.series}, series_CPI);
    CPI = T{:, data(idx).VariableNames};
    Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    Pi3 = [NaN(3, 1);...
        T{4:end,   data(idx).VariableNames}./...
        T{1:end-3, data(idx).VariableNames}];
    Pi12 = [NaN(12, 1);...
        T{13:end,   data(idx).VariableNames}./...
        T{1:end-12, data(idx).VariableNames}];
    
    
    if log_approx

        Pi = log(Pi)*100;
        Pi3 = log(Pi3)*100;
        Pi12 = log(Pi12)*100;

    else

        Pi = (Pi -1)*100;
        Pi3 = (Pi3 -1)*100;
        Pi12 = (Pi12 -1)*100;

    end



    CPI = CPI(idx_sample);
    Pi = Pi(idx_sample);
    Rf = Rf(idx_sample);
    Rf_lag = Rf_lag(idx_sample);
    Rf_lead = Rf_lead(idx_sample);
    
    
    N = size(Pi, 1)- 1;
        
    % estimate equation (6)
    Y = Pi(1:end-1);
    X = [ones(N, 1), Rf_lead(1:end-1)];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    %tab1_ = [N, B', sqrt(diag(Sigma_hat)'), STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:4, 8, 12])];

    Pi_hat = [1, Rf_lead(end)]*B;
    err = Pi(end) - Pi_hat;
    TB_model_results_ = [Pi_hat, err, N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:4, 8, 12])];

    TB_model_results = [TB_model_results; TB_model_results_];    
    %{
    EITB = [T(idx_sample, 'date'),...
        array2table(X*B, 'VariableNames', {'EITB'})];
    UITB = [T(idx_sample, 'date'),...
        array2table(Pi-X*B, 'VariableNames', {'UITB'})];
    
    fname = sprintf('EITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB = fullfile(root, 'EITB', fname);
    writetable(EITB, fpath_EITB)
    
    fname = sprintf('UITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB = fullfile(root, 'UITB', fname);
    writetable(UITB, fpath_UITB)
    %}
    
end



VariableNames = {...
    'Pi_hat',...
    'err',...
    'nobs',...
    'constant',...
    'TB',...
    't_stat_constant',...
    't_stat_TB',...
    'R_sq',...
    'R_sq_adj',...
    's(eta)',...
    'rho_1',...
    'rho_2',...
    'rho_3',...
    'rho_4',...
    'rho_8',...
    'rho_12',...
    };
TB_model_results = array2table(TB_model_results, 'VariableNames', VariableNames);
TB_model_results = [sample_tab, TB_model_results];

%tab1.current_time = repmat({datestr(now)}, size(tab1, 1));
%{
clf
subplot(1,2,1)
hold on
plot(TB_model_results.max_date, TB_model_results.Pi_hat)
plot(TB_model_results.max_date, TB_model_results.Pi_hat+TB_model_results.err)
TB_model_results
sqrt(mean(TB_model_results.err.^2, 'omitnan'))
subplot(1,2,2)
plot(TB_model_results.max_date, TB_model_results.R_sq)
return
%}
%TB_model_results.current_time = repmat({datestr(now)}, size(TB_model_results, 1), 1);


fpath_TB_model_results = fullfile(root, 'TB_model_usa.csv');


writetable(TB_model_results, fpath_TB_model_results)




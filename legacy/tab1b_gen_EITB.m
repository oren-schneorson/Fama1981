%{

Table 1 in Fama (1981) specifies two types of regressions:
a. pp. 549
(6) Pi_t = -b_0 - b_1 dA_t - b_2 dRf_t + b_3 dM_t + eta_t
b. pp. 547
(4) Pi_t = alpha_{t-1} + TB_{t-1} + eta_t

In this script I ran the regression
Pi_t = alpha + beta*TB_{t-1} + eta_t

which is different to Fama's wandering intercept regressions. See
gen_EITB_UITB.m for the correct process, and rerun this...
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
matlab_dir_TASE = fullfile(matlab_dir, 'TASE');

lib_data = fullfile('/media', username, 'D', 'data');
lib_israel = fullfile(lib_data, 'Israel');

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
freq = 'M'; % frequency to be used
log_approx = true;


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

%
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
series_Rf = 'RF_FF';
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

    idx = ismember({data.series}, series_ea);
    A = T{:, data(idx).VariableNames};
    dA = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    dA3 = [NaN(3, 1);...
        T{4:end,   data(idx).VariableNames}./...
        T{1:end-3, data(idx).VariableNames}];
    dA12 = [NaN(12,1);...
        T{13:end,   data(idx).VariableNames}./...
        T{1:end-12, data(idx).VariableNames}];

    dA12_f = [dA12(13:end); NaN(12, 1)];

    
    idx = ismember({data.series}, series_ms);
    M = sum(T{:, [data(idx).VariableNames]}, 2);
    dM = [NaN;...
        M(2:end)./...
        M(1:end-1)];
    dM3 = [NaN(3, 1);...
        M(4:end)./...
        M(1:end-3)];
    dM12 = [NaN(12,1);...
        M(13:end)./...
        M(1:end-12)];
    
    
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

        dA = log(dA)*100;
        dA3 = log(dA3)*100;
        dA12 = log(dA12)*100;
        dA12_f = log(dA12_f)*100;

        dM = log(dM)*100;
        dM3 = log(dM3)*100;
        dM12 = log(dM12)*100;

    else

        Pi = (Pi -1)*100;
        Pi3 = (Pi3 -1)*100;
        Pi12 = (Pi12 -1)*100;

        dA = (dA -1)*100;
        dA3 = (dA3 -1)*100;
        dA12 = (dA12 -1)*100;
        dA12_f = (dA12_f -1)*100;

        dM = (dM -1)*100;
        dM3 = (dM3 -1)*100;
        dM12 = (dM12 -1)*100;

    end



    CPI = CPI(idx_sample);
    Pi = Pi(idx_sample);
    Rf = Rf(idx_sample);
    Rf_lag = Rf_lag(idx_sample);
    Rf_lead = Rf_lead(idx_sample);
    dM = dM(idx_sample);
    dA = dA(idx_sample);
    dA12_f = dA12_f(idx_sample);
    
    
    N = size(Pi, 1);
        
    % estimate equation (6)
    Y = Pi;
    X = [ones(N, 1), Rf_lead];
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
    table1b_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:4, 8, 12])];

    table1b = [table1b; table1b_];    

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

    
end



VariableNames = {...
    'nobs',...
    'constant',...
    'TB_lag',...
    't_stat_constant',...
    't_stat_TB_lag',...
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
table1b = array2table(table1b, 'VariableNames', VariableNames);
table1b = [sample_tab, table1b];
%tab1.current_time = repmat({datestr(now)}, size(tab1, 1));
table1b.current_time = repmat({datestr(now)}, size(table1b, 1), 1);


fpath_table1b = fullfile(root, 'Fama1981', sprintf('table1b.%s.csv', freq));


if exist(fpath_table1b, 'file') == 2 && false
    table1b_ = readtable(fpath_table1b);
    table1b_.Properties.VariableNames = table1b.Properties.VariableNames;
    table1b = [table1b; table1b_];
end
%}
%{
% save each iteration alone

fpath_table1b = fullfile(root, 'Fama1981', 'table1b.csv');
counter = 0;
while exist(fpath_table1b, 'file') == 2
    fpath_tab1 = fullfile(root, 'Fama1981', sprintf('table1b_%d.csv', counter));
    counter = counter + 1;
end
        
%}
%tab1 = sortrows(tab1, {'min_date', 'max_date'});
writetable(table1b, fpath_table1b)




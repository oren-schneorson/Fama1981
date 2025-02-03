%{

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
lag_inflation = 0;

%{
% original dates
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
     datetime(1953,1,31), datetime(2007,7,31), datetime(1953,1,31),...
     datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(2019,12,31), datetime(1977,12,31),...
    datetime(1976,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

tab1a = [];
tab2 = [];
tab3 = [];
tab4 = [];


series_CPI = 'CPIAUCSL';
series_Rf = 'RF_FF';
%series_ms = 'BOGMBASE';
series_ms = {'CURRCIR', 'RESBALNS'};
series_ea = 'INDPRO';
series_C = {'PCENDG', 'PCES'};
series_PD = {'PCEPINDG', 'PCEPIS'};


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
        %T_.date = T_.date + calmonths(lag_inflation);
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

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

T = fillmissing(T, 'previous');




% Sample period
for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);


    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames};

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
    dM = dM(idx_sample);
    dA = dA(idx_sample);
    dA12_f = dA12_f(idx_sample);
    
    
    N = size(Pi, 1);
        
    % estimate equation (6)
    Y = Pi;
    X = [ones(N, 1), dM, dA, dA12_f];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    %tab1a_ = [N, B', sqrt(diag(Sigma_hat)'), STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:4, 8, 12])];
    tab1a_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:4, 8, 12])];


    tab1a = [tab1a; tab1a_];

    
end



VariableNames = {...
    'nobs',...
    'constant',...
    'BG',...
    'DPR',...
    'DPR+12',...
    't_stat_constant',...
    't_stat_BG',...
    't_stat_DPR',...
    't_stat_DPR+12',...
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
tab1a = array2table(tab1a, 'VariableNames', VariableNames);
tab1a = [sample_tab, tab1a];
%tab1a.current_time = repmat({datestr(now)}, size(tab1a, 1));
tab1a.current_time = repmat({datestr(now)}, size(tab1a, 1), 1);


fpath_tab1a = fullfile(root, 'Fama1981', 'tab1a.csv');

if exist(fpath_tab1a, 'file') == 2 && false
    tab1a_ = readtable(fpath_tab1a);
    tab1a_.Properties.VariableNames = tab1a.Properties.VariableNames;
    tab1a = [tab1a; tab1a_];
end
%}
%{
% save each iteration alone

fpath_tab1a = fullfile(root, 'Fama1981', 'tab1a.csv');
counter = 0;
while exist(fpath_tab1a, 'file') == 2
    fpath_tab1a = fullfile(root, 'Fama1981', sprintf('tab1a_%d.csv', counter));
    counter = counter + 1;
end
        
%}
%tab1a = sortrows(tab1a, {'min_date', 'max_date'});
writetable(tab1a, fpath_tab1a)




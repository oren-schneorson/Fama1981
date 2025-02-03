%{

Estimating expected inflation from interest rate models.
Fama, 1975. Short-Term Interest Rates as Predictors of Inflation, AER

pp. 270

(3) r_t = R_t + Delta_t + R_t * Delta

pp.272
(15) E[ r_t | phi^m_{t-1} ] = alpha_0 + gamma R_t

(17) E[ pi_t | phi^m_{t-1} ] = alpha_0 + alpha_1 R_t
alpha_1 = gamma - 1

if E[ r_t | phi^m_{t-1} ] = E[ r_t ] (real rate constant through time)
gamma=0, alpha_1 = -1

Run the regressions:
(19) Delta_t = alpha_0 + alpha_1 R_t + epsilon_t
test: gamma == 0, alpha_1 == -1

(21) Delta_t = alpha_0 + alpha_1 R_t +alpha_2 pi_{t-1} + epsilon_t
test: alpha_2 == 0
no autocorrelation in espilon_t for all lags


data:
R_t:
The one-month nominal rate of interest R_t used in the tests is the return
from the end of month t-1 to the end of month t on the Treasury Bill that
matures closest to the end of month t.
The data are from the quote sheets of Salomon Brothers

CPI:
BLS CPI

Sample:
The tests cover the period from January 1953 through July 1971. Tests for
periods prior to 1953 would be meaningless. First, during World War II and
up to the Trea- sury-Federal Reserve Accord of 1951, interest rates on
Treasury Bills were pegged by the government.

Small note:
the exact expression (3) is used to compute r_t in the empirical work.


TODO:
# generate results for longer maturity bonds

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
lib_ir = fullfile(lib_data, 'ir');
lib_pop = fullfile(lib_data, 'population');
lib_consumption = fullfile(lib_data, 'consumption_usa');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
override_period = false;
lag_inflation = 0;

sample_tab = table(...
    [1:8]',...
    [datetime(1953,1,31), datetime(1953,1,31), datetime(2000,1,31),...
     datetime(1953,1,31), datetime(2007,7,31), datetime(1953,1,31),...
     datetime(1954,1,31), datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(2019,12,31), datetime(1977,12,31),...
    datetime(1976,12,31), datetime(1977,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});


tab1 = [];
tab2 = [];
tab3 = [];
tab4 = [];


series_CPI = 'CPIAUCSL';
series_Rf = 'RF_FF';
series_C = {'PCENDG', 'PCES'};
series_PD = {'PCEPINDG', 'PCEPIS'};


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
        % Boons lag inflation one period
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
    
    
    idx = ismember({data.series}, series_CPI);
    CPI = T{:, data(idx).VariableNames};
    Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    Pi12 = [NaN(12,1);...
        T{13:end,   data(idx).VariableNames}./...
        T{1:end-12, data(idx).VariableNames}];
    
    pi_t = 1./CPI;
    Delta_t = [NaN;...
        pi_t(2:end)./pi_t(1:end-1)-1];
    Delta_t_lag = [NaN; Delta_t(1:end-1)];
    R_t = Rf/100;
    
    if log_approx
        Pi = log(Pi)*100;
        Pi12 = log(Pi12)*100;
    else
        Pi = (Pi -1)*100;
        Pi12 = (Pi12 -1)*100;
    end



    CPI = CPI(idx_sample);
    Pi = Pi(idx_sample);
    Delta_t = Delta_t(idx_sample);
    Delta_t_lag = Delta_t_lag(idx_sample);
    Rf = Rf(idx_sample);
    R_t = R_t(idx_sample);
    r_t = R_t + Delta_t + R_t .* Delta_t;
    
    N = size(Pi, 1);

    tab1_ = [...
        arrayfun(@(lag) autocorr_(Delta_t, lag), 1:12)';...
        1/(N-1)^.5;...
        mean(Delta_t);...
        std(Delta_t);...
        N-1];

    tab2_ = [...
        arrayfun(@(lag) autocorr_(r_t, lag), 1:12)';...
        1/(N-1)^.5;...
        mean(r_t);...
        std(r_t);...
        N-1];
        
    % estimate equation (19)
    Y = Delta_t;
    X = [ones(N, 1), R_t];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    % save regression results
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    
    % save regression results
    % not using hac
    tab3_ = [B', sqrt(diag(Sigma_hat))', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), 1:3)];

    % estimate equation (21)
    Y = Delta_t;
    X = [ones(N, 1), R_t, Delta_t_lag];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    sigma_hat = R'*R/(N-K);
    
    Sigma_hat = inv(D)*sigma_hat;

    % save regression results
    % not using hac
    tab4_ = [B', sqrt(diag(Sigma_hat)'), STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), 1:3)];


    tab1 = [tab1, tab1_];
    tab2 = [tab2, tab2_];  
    tab3 = [tab3; tab3_];
    tab4 = [tab4; tab4_];


    
    
end



VariableNames = {...
    'rho_1',...
    'rho_2',...
    'rho_3',...
    'rho_4',...
    'rho_5',...
    'rho_6',...
    'rho_7',...
    'rho_8',...
    'rho_9',...
    'rho_10',...
    'rho_11',...
    'rho_12',...
    'sigma_rho_1',...
    'mean_Delta_t',...
    'std_Delta_t',...
    'T-1',...
    };
tab1 = array2table(tab1', 'VariableNames', VariableNames);


VariableNames = {...
    'rho_1',...
    'rho_2',...
    'rho_3',...
    'rho_4',...
    'rho_5',...
    'rho_6',...
    'rho_7',...
    'rho_8',...
    'rho_9',...
    'rho_10',...
    'rho_11',...
    'rho_12',...
    'sigma_rho_1',...
    'mean_r_t',...
    'std_r_t',...
    'T-1',...
    };
tab2 = array2table(tab2', 'VariableNames', VariableNames);

VariableNames = {...
    'a_0',...
    'a_1',...
    's(a_0)',...
    's(a_1)',...
    'R_sq',...
    'R_sq_adj',...
    's(e)',...
    'rho_1(e)',...
    'rho_2(e)',...
    'rho_3(e)',...
    };
tab3 = array2table(tab3, 'VariableNames', VariableNames);

VariableNames = {...
    'a_0',...
    'a_1',...
    'a_2',...
    's(a_0)',...
    's(a_1)',...
    's(a_2)',...
    'R_sq',...
    'R_sq_adj',...
    's(e)',...
    'rho_1(e)',...
    'rho_2(e)',...
    'rho_3(e)',...
    };


tab4 = array2table(tab4, 'VariableNames', VariableNames);


tab1 = [sample_tab, tab1];
tab2 = [sample_tab, tab2];
tab3 = [sample_tab, tab3];
tab4 = [sample_tab, tab4];


%rows2vars(tab1)
%rows2vars(tab2)




fpath_tab1 = fullfile(root, 'Fama1975', 'tab1.csv');
fpath_tab2 = fullfile(root, 'Fama1975', 'tab2.csv');
fpath_tab3 = fullfile(root, 'Fama1975', 'tab3.csv');
fpath_tab4 = fullfile(root, 'Fama1975', 'tab4.csv');

writetable(tab1, fpath_tab1)
writetable(tab2, fpath_tab2)
writetable(tab3, fpath_tab3)
writetable(tab4, fpath_tab4)



%{
% as in the paper. I inverted tab1 for easier reading by python
VariableNames = strcat('date_',...
    strrep(cellstr(char(sample_tab.min_date)), '-', ''),...
    '__',...
    strrep(cellstr(char(sample_tab.max_date)), '-', '')...
    );

RowNames = {...
    'rho_1',...
    'rho_2',...
    'rho_3',...
    'rho_4',...
    'rho_5',...
    'rho_6',...
    'rho_7',...
    'rho_8',...
    'rho_9',...
    'rho_10',...
    'rho_11',...
    'rho_12',...
    'sigma_rho_1',...
    'mean_Delta_t',...
    'std_Delta_t',...
    'T-1',...
    };


tab1 = array2table(tab1, 'VariableNames', VariableNames', 'RowNames', RowNames);
%tab1.Properties.VariableNames'
fpath = fullfile(root, 'Fama1975', 'tab1.csv');
writetable(tab1, fpath, 'WriteRowNames', true)


%}

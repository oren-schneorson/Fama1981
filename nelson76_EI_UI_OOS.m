%{

This script replicates desc stats in tables 1-2, Nelson (1976)+
extension to include breakdown to expected and unexpected inflation.

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


lib_EITB = fullfile(root, 'EITB');
lib_UITB = fullfile(root, 'UITB');
lib_crsp = fullfile(lib_data, 'fama_french_factors', 'idx_mth');
lib_ir = fullfile(lib_data, 'ir');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
lag_inflation = 1;
log_approx = true;
%log_approx = false;
freq = 'M';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
to_Quarterly = to_Annual;
vintage = true;
%model = 'AR';
model = 'TB';
%datetime.setDefaultFormats('defaultdate','yyyy-MM-dd')
pcent_mult = 1; % set as 100 for percent

%{
if ~vintage
    lib_ns_vintage = lib_ns;
end
%}


%{
% original dates: Nelson 1976
sample_tab = table(...
    [1:4]',...
    [datetime(1953,6,30), datetime(1953,6,30), datetime(1953,6,30), datetime(1964,1,31)]',...
    [datetime(1971,4,30), datetime(1974,2,28), datetime(1963,12,31), datetime(1974,2,28)]',...
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

tab_EI = table;
tab_UI = table;
tab_EI_UI = table;

VariableNames_EI = {'period', 'min_date', 'max_date',...
    'constant', 'beta_EI',...
    'tstat_constant', 'tstat_beta_EI',...
    'DW', 'R_sq', 'R_sq-adj', 'F', 'F_nc'};
VariableNames_UI = {'period', 'min_date', 'max_date',...
    'constant', 'beta_UI',...
    'tstat_constant', 'tstat_beta_UI',...
    'DW', 'R_sq', 'R_sq-adj', 'F', 'F_nc'};
VariableNames_EI_UI = {'period', 'min_date', 'max_date',...
    'constant', 'beta_EI', 'beta_UI',...
    'tstat_constant', 'tstat_beta_EI', 'tstat_beta_UI',...
    'DW', 'R_sq', 'R_sq-adj', 'F', 'F_nc'};
RowNames = {'Const.', 'beta_EI', 'beta_UI'};


series_CPI = 'CPIAUCSL';
%series_CPI = 'CWUR0000SA0';

series_Rf = 'RF_FF';
series_crsp = 'CRSP_vw';



data = struct(...
    'series', {...
    'CWUR0000SA0',... CPI, not seasonally adjusted, BLS
    'CPIAUCSL',... Price index of all urban consumers
    series_crsp,...
    series_Rf,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_crsp,...
    lib_ir,...
    %'',...
    },...
    'src', {...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    %'',...
    },...
    'VariableNames', cell(1,4)...
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
        %{
        if strcmp(data(i).series, 'CWUR0000SA0')
            cpi_nsa = T_.CWUR0000SA0;
            cpi_nsa(T_.date.Year< 1947) = NaN;
            cpi_nsa(T_.date.Year> 1979) = NaN;
            cpi_sa = sa_adj(cpi_nsa, 12);
            T_.CWUR0000SA0 = cpi_sa;
        end
        %}
        %{
        clf
        hold on
        plot(T_.date, log(cpi_nsa))
        plot(T_.date, log(cpi_sa))
        return
        %}
        
    end

    
    if ~exist('T', 'var')
        T = T_;
        metadata = metadata_;
    else
        if any(ismember(T.Properties.VariableNames, data(i).series))
            continue
        end

        %{
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true, 'Type', 'left');
        %}
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true');
        metadata = [metadata; metadata_];
        
    end
    

end


fpath = fullfile(matlab_dir, 'NRC', 'arima_101_usa.csv')
arima_101 = readtable(fpath);
arima_101 = renamevars(arima_101, 'Pi_hat', 'Pi_hat_101');


if strcmp(freq, 'M')
    %pass
elseif strcmp(freq, 'A')
    T = table2timetable(T);
    T = retime(T, 'yearly', to_Annual);
    T = timetable2table(T);
elseif strcmp(freq, 'Q')
    T = table2timetable(T);
    T = retime(T, 'quarterly', to_Quarterly);
    T = timetable2table(T);
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

%T = fillmissing(T, 'previous');

T = outerjoin(T, arima_101(:, {'date', 'Pi', 'Pi_hat_101', 'epsilon_101'}),...
    'type', 'Left', 'Keys', 'date', 'MergeKeys', true);

fpath = fullfile(root, 'FGLS_OOS.xlsx');
TB_model = readtable(fpath);
%TB_model = renamevars(TB_model, 'Pi_hat', 'Pi_hat_TB');
%TB_model = renamevars(TB_model, 'err', 'epsilon_TB');

TB_model

clf

return
T = outerjoin(T, TB_model(:, {'date', 'Pi_hat_TB', 'epsilon_TB'}),...
    'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
%{
% old version of TB model, in-sample, no wandering intercept correction
fpath = fullfile(root, 'TB_model_usa.csv');
TB_model = readtable(fpath);
TB_model = renamevars(TB_model, 'max_date', 'date');
TB_model = renamevars(TB_model, 'Pi_hat', 'Pi_hat_TB');
TB_model = renamevars(TB_model, 'err', 'epsilon_TB');

T = outerjoin(T, TB_model(:, {'date', 'Pi_hat_TB', 'epsilon_TB'}),...
    'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
%}


% Sample period
%for period = sample_tab.period(end-2)'
for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);

    EIAR = T.Pi_hat_101;
    UIAR = T.epsilon_101;
    

    if strcmp(freq, 'A')
        EITB = EITB*12;
        UITB = UITB*12;
        EIAR = EIAR*12;
        UIAR = UIAR*12;
    elseif strcmp(freq, 'Q')
        EITB = EITB*4;
        UITB = UITB*4;
        EIAR = EIAR*4;
        UIAR = UIAR*4;
    elseif strcmp(freq, 'M')
        % pass
    else
        error('freq must be A, Q or M')
    end
        
    % ***************************************
    % CRSP value weighted inedx

    idx = ismember({data.series}, series_crsp);
    CRSP = T{:, data(idx).VariableNames};
    
    % ***************************************
    % Risk free rate

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames};

    % ***************************************
    % CPI

    idx = ismember({data.series}, series_CPI);
    CPI = T{:, data(idx).VariableNames};
    Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    
    
    dCRSP = [NaN;...
        CRSP(2:end, :)./...
        CRSP(1:end-1, :)];
    
   
    if log_approx
        
        %Pi = log(Pi)*pcent_mult;
        dCRSP = log(dCRSP)*pcent_mult;

    else

        %Pi = (Pi - 1)*pcent_mult;
        dCRSP = (dCRSP- 1)*pcent_mult;

    end
    
    Pi = (Pi - 1)*pcent_mult;
    RS = dCRSP-Pi; % real stock return
    
    EITB = EITB/100*pcent_mult;
    UITB = UITB/100*pcent_mult;
    UIAR = UIAR/100*pcent_mult;
    EIAR = EIAR/100*pcent_mult;
    
    
    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;
    
    if strcmp(model, 'TB')
        EI = EITB;
        UI = UITB;
    elseif strcmp(model, 'AR')
        EI = EIAR;
        UI = UIAR;
    else
        error('model must be AR or TB')
    end
    

    for lead_lag = 1

    Y = dCRSP;
    X = [EI];
    
    if lead_lag < 0
        X = [NaN(abs(lead_lag), size(X,2)); X(1:end-abs(lead_lag),:)];
    else
        X = [X(lead_lag+1:end,:); NaN(lead_lag, size(X,2))];
    end    

    
    X = [ones(size(X, 1), 1), X];
    
    % regressions for EI
    idx = idx_sample & all(~isnan([X, Y]), 2);
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat

    B = [B; NaN];   
    tstat_B = B./[STDB; NaN];
    tstat_B = mat2cell(tstat_B, ones(size(tstat_B)), size(tstat_B,2));
    B_cell = mat2cell(B, ones(size(B)), size(B,2));

    tab_row = table(...
        lead_lag, period, min(T.date(idx_sample)), max(T.date(idx_sample)),...
        B_cell{:}, tstat_B{:}, DW, STATS(1), R_sq_adj, F, STATS(2),...
        'VariableNames', [{'lag'}, VariableNames_EI_UI]);
    
    tab_EI = [tab_EI; tab_row];

    
    % regressions for UI
    Y = dCRSP;
    X = [UI];
    
    if lead_lag < 0
        X = [NaN(abs(lead_lag), size(X,2)); X(1:end-abs(lead_lag),:)];
    else
        X = [X(max(1, lead_lag+1):end,:); NaN(lead_lag, size(X,2))];
    end    
    
    X = [ones(size(X, 1), 1), X];
    
    % regressions for UITB
    idx = idx_sample & all(~isnan([X, Y]), 2);
    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);

    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat

    B = [B(1); NaN; B(2)];
    tstat_B = B./[STDB(1); NaN; STDB(2)];
    
    tstat_B = mat2cell(tstat_B, ones(size(tstat_B)), size(tstat_B,2));
    B_cell = mat2cell(B, ones(size(B)), size(B,2));

    tab_row = table(...
        lead_lag, period, min(T.date(idx_sample)), max(T.date(idx_sample)),...
        B_cell{:}, tstat_B{:}, DW, STATS(1), R_sq_adj, F, STATS(2),...
        'VariableNames', [{'lag'}, VariableNames_EI_UI]);
    
    tab_UI = [tab_UI; tab_row];

    % regressions for EI and UI
    Y = dCRSP;
    X = [EI, UI];
    X = [ones(size(X, 1), 1), X];
    
    if lead_lag < 0
        X = [NaN(abs(lead_lag), size(X,2)); X(1:end-abs(lead_lag),:)];
    else
        X = [X(max(1, lead_lag+1):end,:); NaN(lead_lag, size(X,2))];
    end    

    
    idx = idx_sample & all(~isnan([X, Y]), 2);

    Y = Y(idx, :);
    X = X(idx, :);
    
    N = size(X, 1);
    K = size(X, 2);

    [B,BINT,R,~,STATS] = regress(Y, X);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    sigma_hat = sqrt(R'*R/(N-K));
    D = X'*X; % information matrix
    
    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));

    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat

    tstat_B = B./STDB;
    tstat_B = mat2cell(tstat_B, ones(size(tstat_B)), size(tstat_B,2));
    B_cell = mat2cell(B, ones(size(B)), size(B,2));

    tab_row = table(...
        lead_lag, period, min(T.date(idx_sample)), max(T.date(idx_sample)),...
        B_cell{:}, tstat_B{:}, DW, STATS(1), R_sq_adj, F, STATS(2),...
        'VariableNames', [{'lag'}, VariableNames_EI_UI]);
    
    tab_EI_UI = [tab_EI_UI; tab_row];

    
    
    
    end    
    
    
end


PycharmProjects_dir = fullfile('/home', username, 'PycharmProjects');

tab_EI_UI.R_sq
fname = sprintf('data_tab_EI%s.%s.csv', model, freq);
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', fname);
writetable(tab_EI, fpath)

fname = sprintf('data_tab_UI%s.%s.csv', model, freq);
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', fname);
writetable(tab_UI, fpath)

fname = sprintf('data_tab_EI%s_UI%s.%s.csv', model, model, freq);
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', fname);
writetable(tab_EI_UI, fpath)


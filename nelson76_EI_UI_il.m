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
lib_data = fullfile('/media', username, 'D', 'data', 'BOI');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')

lib_israel = fullfile(lib_data, 'Israel');
lib_cpi = lib_israel;


lib_EITB = fullfile(root, 'EITB_il', series_Rf);
lib_UITB = fullfile(root, 'UITB_il');
lib_EIPF = lib_israel;
lib_UIPF = lib_israel;

lib_tase_all_share = fullfile(lib_data, 'stock_markets');
%lib_tase_125 = fullfile(lib_israel, 'SM', 'M', 'with_metadata');
lib_tase_125 = fullfile(lib_israel, 'SM_15', 'with_metadata');

lib_ir = fullfile(lib_israel, 'TELBOR', 'M', 'with_metadata');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
lag_inflation = 0;
log_approx = true;
%log_approx = false;
freq = 'M';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
to_Quarterly = to_Annual;
vintage = true;
model = 'PF';
%model = 'AR';
%model = 'TB';
%datetime.setDefaultFormats('defaultdate','yyyy-MM-dd')
pcent_mult = 1; % set as 100 for percent

%{
if ~vintage
    lib_ns_vintage = lib_ns;
end
%}


fpath_EIPF = fullfile(lib_EIPF, 'EIPF.xlsx');
fpath_UIPF = fullfile(lib_EIPF, 'UIPF.xlsx');
fpath_arima = fullfile(matlab_dir, 'NRC', 'arima_101_il.csv')


%{
% for EITB, or EIAR
sample_tab = table(...
    [1:6]',...
    [...
    datetime(1996,3,31), datetime(1996,3,31), datetime(2008,10,31),...
    datetime(2000,1,31), datetime(2000,1,31), datetime(2008,10,31),...
     ]',...
    [...
    datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});

%}

% for EIPF + EIAR
sample_tab = table(...
    [1:9]',...
    [...
    datetime(2001,4,30), datetime(2001,4,30), datetime(2008,10,31),...
    datetime(1996,3,31), datetime(1996,3,31), datetime(2008,10,31),... % EIAR only
    datetime(2001,4,30), datetime(2001,4,30), datetime(2008,10,31),... % EIPF
     ]',...
    [...
    datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    datetime(2024,1,31), datetime(2008,9,30), datetime(2024,1,31),...
    datetime(2024,1,31), datetime(2008,9,30), datetime(2024,1,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});




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


series_CPI = 'CP_NSA.M';

%series_Rf = 'TSB_BAGR_MAKAM_01M.M';
series_Rf = 'BL.TELBOR_01M.M';
series_tase_all_share = 'SPASTT01ILM661N'; % tase-all-share
%series_tase_125 = 'MD_M137.M';
series_tase_125 = 'MD_M137.M_15';
series_sm = series_tase_125;


data = struct(...
    'series', {...
    'CP_NSA.M',... CPI, not seasonally adjusted
    'CP_SA.M',... CPI, seasonally adjusted official series (1996-...)
    series_tase_all_share,...
    series_tase_125,...
    series_Rf,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_tase_all_share,...
    lib_tase_125,...
    lib_ir,...
    %'',...
    },...
    'src', {...
    'FAME',...
    'FAME',...
    'FRED',...
    'FAME',...
    'FAME',...
    %'',...
    },...
    'VariableNames', cell(1,5)...
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

    % robust date, start of month, always first of the month
    aux_date = T_.date+1;
    if strcmp(data(i).series, series_tase_125) && contains(series_tase_125, '_15')
        % from post 15th of the month -->
        % to end of pervious month, to align with unexpected inflation
        T_.date.Day = 1;
        T_.date = T_.date - caldays(1);

    else
        % drop observations that are not either end of month or start of month
        idx = aux_date.Day == 1 | T_.date.Day == 1;
        T_ = T_(idx, :);
    end

    aux_date = T_.date+1;
    if all(aux_date.Day == 1)
        %T_.date = aux_date-calmonths(1);
        T_.date = aux_date;
    end

    if strcmp(data(i).series, series_CPI)
        % Boons lag inflation one period, but Fama does not.
        T_.date = T_.date - calmonths(lag_inflation);
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
        %metadata = [metadata; metadata_];
        
    end
    

end


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

CP_SA = sa_adj(T.CP_NSA0x2EM(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')+12), 12);
T.CP_SA = T.CP_SA0x2EM;
T.CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1) = ...
    CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1);

%{
% check sa_adj valid, pre 1996 data.
clf
hold on
plot(T.date, T.CP_SA)
plot(T.date, T.CP_SA0x2EM)
%plot(T.date, T.CP_NSA0x2EM)
%}

arima_101 = readtable(fpath_arima);

EIAR = arima_101(:, {'date', 'EIAR1_101'});
UIAR = arima_101(:, {'date', 'UIAR1_101'});

EIAR.Properties.VariableNames = {'date', 'EIAR'};
UIAR.Properties.VariableNames = {'date', 'UIAR'};


EIPF = readtable(fpath_EIPF);
UIPF = readtable(fpath_UIPF);

EIPF.date = datetime(EIPF.date);
UIPF.date = datetime(UIPF.date);

%{
% shift dates of EIAR, EIAR one month backward, because EI is t-1.
EIAR.date = EIAR.date + 1;
EIAR.date = EIAR.date - calmonths(1);
EIAR.date = EIAR.date -1;

EIPF.date = EIPF.date + 1;
EIPF.date = EIPF.date - calmonths(1);
EIPF.date = EIPF.date -1;
%}

%{
T = outerjoin(T, arima_101(:, {'date', 'Pi', 'Pi_hat_101', 'epsilon_101'}),...
    'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
%}

%{
% old version of TB model, in-sample, no wandering intercept correction
fpath = fullfile(root, 'TB_model_usa.csv')
TB_model = readtable(fpath);
TB_model = renamevars(TB_model, 'max_date', 'date');
TB_model = renamevars(TB_model, 'Pi_hat', 'Pi_hat_TB');
TB_model = renamevars(TB_model, 'err', 'epsilon_TB');

T = outerjoin(T, TB_model(:, {'date', 'Pi_hat_TB', 'epsilon_TB'}),...
    'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
%}


T = outerjoin(T, EIAR(:, {'date', 'EIAR'}),...
    'Keys', {'date'}, 'MergeKeys', true, 'Type', 'left');
T = outerjoin(T, UIAR(:, {'date', 'UIAR'}),...
    'Keys', {'date'}, 'MergeKeys', true, 'Type', 'left');
T = outerjoin(T, EIPF(:, {'date', 'EIPF'}),...
    'Keys', {'date'}, 'MergeKeys', true, 'Type', 'left');
T = outerjoin(T, UIPF(:, {'date', 'UIPF'}),...
    'Keys', {'date'}, 'MergeKeys', true, 'Type', 'left');


% Sample period
%for period = sample_tab.period(end-2)'
%for period = sample_tab.period(4:end)'
for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);
    %{
    series_EITB = sprintf('EITB__%s-%s.%s',...
        min_date, max_date, freq);
    series_UITB = sprintf('UITB__%s-%s.%s',...
        min_date, max_date, freq);
    
    fname = [series_EITB, '.csv'];
    fpath_EITB = fullfile(lib_EITB, fname);
    EITB = readtable(fpath_EITB);

    fname = [series_UITB, '.csv'];
    fpath_UITB = fullfile(lib_UITB, fname);
    UITB = readtable(fpath_UITB);

    T = outerjoin(T, EITB,...
        'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
    T = outerjoin(T, UITB,...
        'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
    
    
    EITB = T.EITB;
    UITB = T.UITB;
    T = removevars(T, {'EITB', 'UITB'});  
    %}



    % old version of TB model, no wandering intercept
    %EITB_period = T.Pi_hat_TB;
    %UITB_period = T.epsilon_TB;
    EIPF_period = T.EIPF;
    UIPF_period = T.UIPF;
    EIAR_period = T.EIAR;
    UIAR_period = T.UIAR;
    

    if strcmp(freq, 'A')
        %EITB_period = EITB_period*12;
        %UITB_period = UITB_period*12;
        EIPF_period = EIPF_period*12;
        UIPF_period = UIPF_period*12;
        EIAR_period = EIAR_period*12;
        UIAR_period = UIAR_period*12;
    elseif strcmp(freq, 'Q')
        %EITB_period = EITB_period*4;
        %UITB_period = UITB_period*4;
        EIPF_period = EIPF_period*4;
        UIPF = UIPF_period*4;
        EIAR_period = EIAR_period*4;
        UIAR_period = UIAR_period*4;
    elseif strcmp(freq, 'M')
        % pass
    else
        error('freq must be A, Q or M')
    end
        
    % ***************************************
    % TASE value weighted inedx

    idx = ismember({data.series}, series_sm);
    TASE = T{:, data(idx).VariableNames};
    
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
    
    
    dTASE = [NaN;...
        TASE(2:end, :)./...
        TASE(1:end-1, :)];
    
   
    if log_approx
        
        %Pi = log(Pi)*pcent_mult;
        dTASE = log(dTASE)*pcent_mult;

    else

        %Pi = (Pi - 1)*pcent_mult;
        dTASE = (dTASE- 1)*pcent_mult;

    end
    
    Pi = (Pi - 1)*pcent_mult;
    RS = dTASE-Pi; % real stock return
    
    %EITB_period = EITB_period/100*pcent_mult;
    %UITB_period = UITB_period/100*pcent_mult;
    EIPF_period = EIPF_period/100*pcent_mult;
    UIPF_period = UIPF_period/100*pcent_mult;
    UIAR_period = UIAR_period/100*pcent_mult;
    EIAR_period = EIAR_period/100*pcent_mult;
    
    
    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;
    
    if strcmp(model, 'TB')
        EI = EITB_period;
        UI = UITB_period;
    elseif strcmp(model, 'AR')
        EI = EIAR_period;
        UI = UIAR_period;
    elseif strcmp(model, 'PF')
        EI = EIPF_period;
        UI = UIPF_period;
    else
        error('model must be AR, TB or PF')
    end
    

    for lead_lag = 1

    Y = dTASE;
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
    Y = dTASE;
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
    Y = dTASE;
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



tab_EI_UI



PycharmProjects_dir = fullfile('/home', username, 'PycharmProjects');

fname = sprintf('data_tab_EI%s.%s.csv', model, freq);
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', 'il', fname);
writetable(tab_EI, fpath)

fname = sprintf('data_tab_UI%s.%s.csv', model, freq);
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', 'il', fname);
writetable(tab_UI, fpath)

fname = sprintf('data_tab_EI%s_UI%s.%s.csv', model, model, freq);
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', 'il', fname);
writetable(tab_EI_UI, fpath)


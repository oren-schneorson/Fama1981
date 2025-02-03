%{

This plots the CX/NS and ROC series for Israel

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

lib_cpi = fullfile(lib_israel);
lib_ms = fullfile(lib_israel, 'monagg');
lib_mbase = fullfile(lib_israel, 'money_base');
lib_reserves = fullfile(lib_israel, 'reserves');
lib_roc = fullfile(lib_israel);


lib_gnp = fullfile(lib_israel, 'GNP');  % gross national product
lib_capital_stock = fullfile(lib_israel, 'capital_stock');
lib_cx = fullfile(lib_capital_stock, 'inv');  % capital expenditure (investment)
lib_ns = fullfile(lib_capital_stock, 'nstock');  % net capital stock

%lib_ea = fullfile(lib_israel, 'indprod');  % economic activity
lib_ea = fullfile(lib_israel);  % economic activity, monthly series

%lib_ir = fullfile(lib_israel);
lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M', 'with_metadata');
lib_pop = fullfile(lib_israel, 'POP');
lib_consumption = fullfile(lib_israel, 'C');
%lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');

lib_EITB = fullfile(root, 'EITB_il');
lib_UITB = fullfile(root, 'UITB_il');



% flag to override period in downstream scripts: false by default
lag_inflation = 0;
%log_approx = true;
log_approx = false;
freq = 'Q';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
pcent_mult = 1; % set as 100 for percent


sample_tab = table(...
    [1:6]',...
    [...
    datetime(1999,12,31), datetime(1999,12,31), datetime(2008,10,31),...
    datetime(2001,12,31), datetime(2001,12,31), datetime(2008,10,31),...
     ]',...
    [...
    datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});

sample_tab = sample_tab(4:end, :);

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});


EI_measure = {'AR'};
series_CPI = 'CP_SA.M';
series_CPI_NSA = 'CP_NSA.M';

if strcmp(freq, 'M')
    series_Rf = 'TSB_BAGR_MAKAM_01M.M'; % 1 month MAKAM, annualized
elseif strcmp(freq, 'Q')
    series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized
elseif strcmp(freq, 'A')
    series_Rf = 'TSB_BAGR_MAKAM_12M.M'; % 12 month MAKAM, annualized
else
    error('freq must be M, Q or A')
end

series_ms = 'ZA108.M'; % waiting for Elad from MOS
%series_ms = {'A346.M_E', 'ZA215.M'}; % ask Elad about RESNALNS, 

series_ea = 'mdd13.level.m'; % excl diamonds
%series_ea = 'CLS11.TPR_C.Q_SA_CHAINED'; % excl diamonds

series_tax = 'BOGZ1FL143178005Q';
series_mip = 'BOGZ1FA106130001Q'; % monetary interest paid
series_pbt = 'A464RC1Q027SBEA'; % profit before tax
series_pbt_tr = 'BOGZ1FU106060005Q'; % profit before tax
series_pbt_tr_incl = 'BOGZ1FU106060035Q'; % profit before tax, incl iva and ccadj
series_iva = 'NCBIVDQ027S';

series_gdp_bs = 'WRK.GDP_BS.Q_N';
series_comp_bs = 'COMP_BS.Q';

% note: there are discrepancies between NIPA and FRED. use NIPA
series_cx = 'AM_INV_TOTALT_Q_N';

%series_gnp = 'GNP.Q_N';
series_gnp = 'GNP.Q_FP';

% residential and nonresidential equipment and structures, capital stock
series_ns = {'AM_NSTOCK_TOTALT_Q_N'}; % includes inventories... not exactly like Fama
series_roc = 'ROC_BS.Q';
series_roc_denom = {'AM_NSTOCK_TOTALT_Q_N'}; % series_ns_fred + inventories
series_C = {'C_@DUR.Q_FP_SA'}; % excl durable, should contain services
% not needed
%series_PD = {'PCEPINDG', 'PCEPIS'}; % I have the real series, don't really
%need this. (implicit deflator)
%series_pop = 'DEM_POP_AVE_CHN.Q';


data = struct(...
    'series', {...
    series_CPI,... General Consumer Price Index, seasonally adjusted
    series_CPI_NSA,... General Consumer Price Index, non seasonally adjusted
    series_Rf,...
    'ZA108.M',...
    'A346.M_E',...
    'ZA215.M',...
    series_gnp,... % real GNP
    series_ea,...
    series_cx,...
    'AM_NSTOCK_ISKIM_Q_N',... % Business Sector; Equipment
    'AM_NSTOCK_ISKIS_Q_N',... % Business Sector; Structures
    'AM_NSTOCK_BILDRES_Q_N',...  % Construction Sector; Residential All
    'AM_NSTOCK_TOTALT_Q_N',...
    series_gdp_bs,... % business sector GDP
    series_comp_bs,... % business sector compensation
    series_roc,...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_ir,...
    lib_mbase,...
    lib_mbase,...
    lib_reserves,...
    lib_gnp,...
    lib_ea,...
    lib_cx,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_roc,...
    lib_roc,...
    lib_roc,...
    %'',...
    },...
    'src', {...
    'FAME',...
    'FAME',...
    'FAME',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'FAME',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'FAME',...
    'FAME',...
    'FAME',...
    },...
    'VariableNames', cell(1,16)...
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
        % T_.date = T_.date + calmonths(lag_inflation);
    end
    
    if ~exist('T', 'var')
        T = T_;
        metadata = metadata_;
    else
        %{
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true, 'Type', 'left');
        %}
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true');
        metadata = [metadata; metadata_];
    end
    

end


% In-sample UIAR
fpath_factor = fullfile(matlab_dir_TASE, 'arima_101_il.csv');
arima_101 = readtable(fpath_factor);
fpath_factor = fullfile(matlab_dir_TASE, 'arima_011_il.csv');
arima_011 = readtable(fpath_factor);

% Out-of-sample UITB
%fpath = '../Fama1981/UITB_il/UITB__1996-03-31-2019-12-31.M.csv';
fpath_factor = fullfile(lib_israel, 'UITB.xlsx');
UITB = readtable(fpath_factor);
%fpath = '../Fama1981/UITB_il/UITB1__1996-03-31-2023-08-31.M.csv';
fpath_factor = fullfile(lib_israel, 'UITB1.xlsx');
UITB1 = readtable(fpath_factor);

fpath_factor = fullfile(lib_israel, 'UIPF_1.0.xlsx');
UIPF = readtable(fpath_factor);
UIPF.date = datetime(UIPF.date);


EITB = renamevars(UITB, 'UITB', 'EITB');
EITB1 = renamevars(UITB1, 'UITB1', 'EITB1');
EIPF = renamevars(UIPF, 'UIPF', 'EIPF');
UIAR = arima_101(:, {'date', 'UIAR1_101'});
EIAR = arima_101(:, {'date', 'EIAR1_101'});

EITB.EITB = UITB.Pi-UITB.UITB;
EIPF.EIPF = UIPF.Pi_SA-UIPF.UIPF;
EITB1.EITB1 = UITB1.Pi-UITB1.UITB1;


% shift expectation one period backward
EITB.date = EITB.date + 1;
EITB.date = EITB.date - calmonths(1);
EITB.date = EITB.date - 1;

EIPF.date = EIPF.date + 1;
EIPF.date = EIPF.date - calmonths(1);
EIPF.date = EIPF.date - 1;

EITB1.date = EITB1.date + 1;
EITB1.date = EITB1.date - calmonths(1);
EITB1.date = EITB1.date - 1;

EIAR.date = EIAR.date+1;
EIAR.date = EIAR.date - calmonths(1);
EIAR.date = EIAR.date-1;



if strcmp(freq, 'M')
    %pass
elseif strcmp(freq, 'Q')
    T = table2timetable(T);
    T = retime(T, 'quarterly', 'lastvalue');
    T = timetable2table(T);

    EITB = table2timetable(EITB);
    EITB = retime(EITB, 'quarterly', 'mean');
    EITB = timetable2table(EITB);
    EITB.date = EITB.date - 1;
    EITB.date = EITB.date + calmonths(3);

    EITB1 = table2timetable(EITB1);
    EITB1 = retime(EITB1, 'quarterly', 'mean');
    EITB1 = timetable2table(EITB1);
    EITB1.date = EITB1.date - 1;
    EITB1.date = EITB1.date + calmonths(3);

    EIAR = table2timetable(EIAR);
    EIAR = retime(EIAR, 'quarterly', 'mean');
    EIAR = timetable2table(EIAR);
    EIAR.date = EIAR.date - 1;
    EIAR.date = EIAR.date + calmonths(3);

    EIPF = table2timetable(EIPF);
    EIPF = retime(EIPF, 'quarterly', 'mean');
    EIPF = timetable2table(EIPF);
    EIPF.date = EIPF.date - 1;
    EIPF.date = EIPF.date + calmonths(3);

    % Unexpected inflation
    UITB = table2timetable(UITB);
    UITB = retime(UITB, 'quarterly', 'mean');
    UITB = timetable2table(UITB);
    UITB.date = UITB.date - 1;
    UITB.date = UITB.date + calmonths(3);

    UITB1 = table2timetable(UITB1);
    UITB1 = retime(UITB1, 'quarterly', 'mean');
    UITB1 = timetable2table(UITB1);
    UITB1.date = UITB1.date - 1;
    UITB1.date = UITB1.date + calmonths(3);

    UIAR = table2timetable(UIAR);
    UIAR = retime(UIAR, 'quarterly', 'mean');
    UIAR = timetable2table(UIAR);
    UIAR.date = UIAR.date - 1;
    UIAR.date = UIAR.date + calmonths(3);

    UIPF = table2timetable(UIPF);
    UIPF = retime(UIPF, 'quarterly', 'mean');
    UIPF = timetable2table(UIPF);
    UIPF.date = UIPF.date - 1;
    UIPF.date = UIPF.date + calmonths(3);


elseif strcmp(freq, 'A')
    T = table2timetable(T);
    T = retime(T, 'yearly', 'lastvalue');
    T = timetable2table(T);

    EITB = table2timetable(EITB);
    EITB = retime(EITB, 'yearly', 'mean');
    EITB = timetable2table(EITB);
    EITB.date = EITB.date - 1;
    EITB.date = EITB.date + calmonths(12);

    EITB1 = table2timetable(EITB1);
    EITB1 = retime(EITB1, 'yearly', 'mean');
    EITB1 = timetable2table(EITB1);
    EITB1.date = EITB1.date - 1;
    EITB1.date = EITB1.date + calmonths(12);

    EIAR = table2timetable(EIAR);
    EIAR = retime(EIAR, 'yearly', 'mean');
    EIAR = timetable2table(EIAR);
    EIAR.date = EIAR.date - 1;
    EIAR.date = EIAR.date + calmonths(12);

    EIPF = table2timetable(EIPF);
    EIPF = retime(EIPF, 'yearly', 'mean');
    EIPF = timetable2table(EIPF);
    EIPF.date = EIPF.date - 1;
    EIPF.date = EIPF.date + calmonths(12);

    % Unexpected inflation
    UITB = table2timetable(UITB);
    UITB = retime(UITB, 'yearly', 'mean');
    UITB = timetable2table(UITB);
    UITB.date = UITB.date - 1;
    UITB.date = UITB.date + calmonths(12);

    UITB1 = table2timetable(UITB1);
    UITB1 = retime(UITB1, 'yearly', 'mean');
    UITB1 = timetable2table(UITB1);
    UITB1.date = UITB1.date - 1;
    UITB1.date = UITB1.date + calmonths(12);

    UIAR = table2timetable(UIAR);
    UIAR = retime(UIAR, 'yearly', 'mean');
    UIAR = timetable2table(UIAR);
    UIAR.date = UIAR.date - 1;
    UIAR.date = UIAR.date + calmonths(12);

    UIPF = table2timetable(UIPF);
    UIPF = retime(UIPF, 'yearly', 'mean');
    UIPF = timetable2table(UIPF);
    UIPF.date = UIPF.date - 1;
    UIPF.date = UIPF.date + calmonths(12);



else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

T = fillmissing(T, 'previous');


% Sample period
for period = sample_tab.period'

    min_date = sample_tab.min_date(period==sample_tab.period);
    max_date = sample_tab.max_date(period==sample_tab.period);

    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    
    %{
    series_EITB = sprintf('EITB__%s-%s.M',...
        min_date, max_date);
    series_UITB = sprintf('UITB__%s-%s.M',...
        min_date, max_date);
    
    fname = [series_EITB, '.csv'];
    fpath_EITB = fullfile(lib_EITB, fname);
    EITB = readtable(fpath_EITB);

    fname = [series_UITB, '.csv'];
    fpath_UITB = fullfile(lib_UITB, fname);
    UITB = readtable(fpath_UITB);
    %}
    T = outerjoin(T, EIAR, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, UIAR, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, EIPF(:, {'date', 'EIPF'}), 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, UIPF(:, {'date', 'UIPF'}), 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, EITB, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, UITB, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    
    EIPF_ = T.EIPF;
    UIPF_ = T.UIPF;
    EIAR_ = T.EIAR1_101;
    UIAR_ = T.UIAR1_101;
    EITB_ = T.EITB;
    UITB_ = T.UITB;
    T = removevars(T, {'EITB', 'UITB', 'EIPF', 'UIPF', 'EIAR1_101', 'UIAR1_101'});  

    eval(sprintf('EI = EI%s_;', char(EI_measure)))
    eval(sprintf('UI = UI%s_;', char(EI_measure)))


    if strcmp(freq, 'A')
        EI = EI*12;
        UI = UI*12;
    elseif strcmp(freq, 'Q')
        EI = EI*4;
        UI = UI*4;
    elseif strcmp(freq, 'M')
        % pass
    else
        error('freq must be A, Q or M')
    end
    % ***************************************
    % Capital expenditures

    idx = ismember({data.series}, series_cx);
    CX = sum(T{:, [data(idx).VariableNames]}, 2, 'omitnan');
    
    % ***************************************
    % Net capital stock

    idx = ismember({data.series}, series_ns);
    NS = sum(T{:, [data(idx).VariableNames]}, 2);

    % ***************************************
    % Return on Capital, ROC
    
    % not seasonally adjusted
    idx = ismember({data.series}, series_roc);
    ROC = sum(T{:, [data(idx).VariableNames]}, 2);

    idx = ismember({data.series}, series_comp_bs);
    COMP_BS = sum(T{:, [data(idx).VariableNames]}, 2);

    idx = ismember({data.series}, series_gdp_bs);
    GDP_BS = sum(T{:, [data(idx).VariableNames]}, 2);

    idx_nan_COMP_BS = isnan(COMP_BS);
    idx_nan_GDP_BS = isnan(GDP_BS);
    idx_nan_ROC = isnan(ROC);
    if strcmp(freq, 'M')
        COMP_BS = [NaN(sum(idx_nan_COMP_BS), 1); sa_adj(COMP_BS(~idx_nan_COMP_BS), 12)];
        GDP_BS = [NaN(sum(idx_nan_GDP_BS), 1); sa_adj(GDP_BS(~idx_nan_GDP_BS), 12)];
        ROC = [NaN(sum(idx_nan_ROC), 1); sa_adj(ROC(~idx_nan_ROC), 12)];
    elseif strcmp(freq, 'Q')
        COMP_BS = [NaN(sum(idx_nan_COMP_BS), 1); sa_adj(COMP_BS(~idx_nan_COMP_BS), 4)];
        GDP_BS = [NaN(sum(idx_nan_GDP_BS), 1); sa_adj(GDP_BS(~idx_nan_GDP_BS), 4)];
        %ROC = [NaN(sum(idx_nan_ROC), 1); sa_adj(ROC(~idx_nan_ROC), 4)];
    end
    %{
    ROC_ = 1-COMP_BS./GDP_BS;
    clf
    hold on
    plot(T.date(idx_sample), ROC(idx_sample))
    plot(T.date(idx_sample), ROC_(idx_sample))
    return
    %}


    % ***************************************
    % Risk free rate

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames};

    % ***************************************
    % Real GNP
    
    idx = ismember({data.series}, series_gnp);
    RGNP = T{:, data(idx).VariableNames};
    dRGNP = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    
    % ***************************************
    % Industrial production
    
    idx = ismember({data.series}, series_ea);
    PR = T{:, data(idx).VariableNames};
    dPR = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];

    % ***************************************
    % Money base
    
    idx = ismember({data.series}, series_ms);
    M = sum(T{:, [data(idx).VariableNames]}, 2);
    dM = [NaN;...
        M(2:end)./...
        M(1:end-1)];
    
    % ***************************************
    % CPI
    if strcmp(freq, 'A')
        idx = ismember({data.series}, series_CPI_NSA);
    else
        idx = ismember({data.series}, series_CPI);
    end

    CPI = T{:, data(idx).VariableNames};
    Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    
        
    CX_NS = CX ./ NS;
    
    %{
    dROC = [NaN;...
        ROC(2:end)./...
        ROC(1:end-1)];
    
    
    dCX_NS = [NaN;...
        CX_NS(2:end)./...
        CX_NS(1:end-1)];    
    %}
    dROC = [NaN;...
        ROC(2:end)-...
        ROC(1:end-1)];
    
    
    dCX_NS = [NaN;...
        CX_NS(2:end)-...
        CX_NS(1:end-1)];    
    

   
    if log_approx
        Pi = log(Pi)*pcent_mult;
        dPR = log(dPR)*pcent_mult;
        dRGNP = log(dRGNP)*pcent_mult;
        dM = log(dM)*pcent_mult;
        %dCX_NS = log(dCX_NS)*pcent_mult;
        %dROC = log(dROC)*pcent_mult;

    else
        Pi = (Pi - 1)*pcent_mult;
        dPR = (dPR - 1)*pcent_mult;
        dRGNP = (dRGNP- 1)*pcent_mult;
        dM = (dM - 1)*pcent_mult;
        %dCX_NS = (dCX_NS - 1)*pcent_mult;
        %dROC = (dROC- 1)*pcent_mult;

    end
    
    EI = EI/100*pcent_mult;
    UI = UI/100*pcent_mult;

    % define leads and lags
    % T = sortrows(T, 'date', 'ascending');  % make sure T is ascending
    dPR_lag = [NaN; dPR(1:end-1)];
    dM_lag = [NaN; dM(1:end-1)];
    dCX_NS_lag = [NaN; dCX_NS(1:end-1)];
    dROC_lag = [NaN; dROC(1:end-1)];
    dRGNP_lag = [NaN; dRGNP(1:end-1)];
    Rf_lag = [NaN; Rf(1:end-1)];
    EI_lag = [NaN; EI(1:end-1)];
    UI_lag = [NaN; UI(1:end-1)];

    dPR_lead = [dPR(2:end); NaN];
    dM_lead = [dM(2:end); NaN];
    dCX_NS_lead = [dCX_NS(2:end); NaN];
    dROC_lead = [dROC(2:end); NaN];
    dRGNP_lead = [dRGNP(2:end); NaN];
    Rf_lead = [Rf(2:end); NaN];
    

    CPI = CPI(idx_sample);
    Pi = Pi(idx_sample);
    Rf = Rf(idx_sample);
    Rf_lag = Rf_lag(idx_sample);
    
    dM = dM(idx_sample);
    
    dPR = dPR(idx_sample);
    dPR_lead = dPR_lead(idx_sample);
    dPR_lag = dPR_lag(idx_sample);

    dRGNP = dRGNP(idx_sample);
    dRGNP_lead = dRGNP_lead(idx_sample);
    dRGNP_lag = dRGNP_lag(idx_sample);
    
    CX = CX(idx_sample);
    NS = NS(idx_sample);
    CX_NS = CX_NS(idx_sample);
    dCX_NS = dCX_NS(idx_sample);
    dCX_NS_lead = dCX_NS_lead(idx_sample);
    dCX_NS_lag = dCX_NS_lag(idx_sample);

    ROC = ROC(idx_sample);
    dROC = dROC(idx_sample);
    dROC_lead = dROC_lead(idx_sample);
    dROC_lag = dROC_lag(idx_sample);

    EI = EI(idx_sample);
    EI_lag = EI_lag(idx_sample);
    UI = UI(idx_sample);
    UI_lag = UI_lag(idx_sample);

    dates = T.date(idx_sample);

    %
    if period == 4
        % CX_NS_il.png
        clf
        hold on
        plot(dates, log(CX), 'b')
        plot(dates, log(NS), 'r')
        yyaxis right
        plot(dates, CX_NS, 'k--')

        legend({'log(CX)', 'log(NS)', 'CX/NS (right)'}, 'Location', 'northwest')
        set(gcf, 'Color', 'w')
        ax = gca;
        ax.YAxis(2).Color = 'k';
        export_fig(fullfile(root, 'gfx', 'CX_NS_il.png'))

        % CX_NS_CI_il.png
        clf
        hold on
        plot(dates, CX_NS, 'k--')
        yyaxis right
        plot(dates, dPR, 'g-')

        legend({'CX/NS', 'Composite Index (CI)'}, 'Location', 'northwest')
        set(gcf, 'Color', 'w')
        ax = gca;
        ax.YAxis(2).Color = 'k';
        export_fig(fullfile(root, 'gfx', 'CX_NS_CI_il.png'))

        clf
        GDP_BS_sa = sa_adj(GDP_BS(idx_sample), 4);
        COMP_BS_sa = sa_adj(COMP_BS(idx_sample), 4);

        hold on
        plot(dates, log(COMP_BS(idx_sample)), 'Color', [0, 0, .5])
        plot(dates, log(COMP_BS_sa), 'b')
        plot(dates, log(GDP_BS(idx_sample)), 'Color', [.7, 0, 0])
        plot(dates, log(GDP_BS_sa), 'r')
        
        yyaxis right
        plot(dates, ROC, 'LineStyle', '-', 'Color', [0.5, 0.5, 0.5])
        plot(dates, 1-COMP_BS_sa./GDP_BS_sa, 'k--')

        legend(str4fig({'log(COMP_BS)', 'log(COMP_BS) (SA)', 'log(GDP_BS)', 'log(GDP_BS) (SA)', '1-ROC', '1-ROC (SA)'}), 'Location', 'northwest')
        set(gcf, 'Color', 'w')
        export_fig(fullfile(root, 'gfx', 'ROC_il.png'))
        ax = gca;
        ax.YAxis(2).Color = 'k';
        export_fig(fullfile(root, 'gfx', 'ROC_il.png'))

        return

        % ROC_CI_il.png
        clf
        hold on
        plot(dates, dROC, 'k--')
        yyaxis right
        plot(dates, dPR, 'g-')

        legend({'dROC (SA)', 'Composite Index (CI)'}, 'Location', 'northwest')
        set(gcf, 'Color', 'w')
        ax = gca;
        ax.YAxis(2).Color = 'k';
        export_fig(fullfile(root, 'gfx', 'ROC_CI_il.png'))

        return
    end  


    
end


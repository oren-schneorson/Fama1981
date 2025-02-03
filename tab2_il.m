%{
Data for Israel.
# No breakdown to nonfinancial/financial)
# Can see Business sector, but not corporate (private) sector, incl.
govermental agencies.


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
lag_inflation = 1;
log_approx = false;
pcent_mult = 1; % set as 100 if want percent
freq = 'A'; % frequency to be used



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

%{
sample_tab = table(...
    [1:3]',...
    [datetime(1996,3,31), datetime(1996,2,31), datetime(2008,10,31),...
     ]',...
    [datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%{
% no EITB for pre 1992, because no short rate
sample_tab = table(...
    [1:3]',...
    [datetime(1980,1,31), datetime(1980,3,31), datetime(2009,1,31),...
     ]',...
    [datetime(2019,12,31), datetime(2008,12,31), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

table2 = [];

series_CPI = 'CP_SA.M';
series_CPI_NSA = 'CP_NSA.M';
EI_measure = {'AR'};

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

% note: there are discrepancies between NIPA and FRED. use NIPA
series_cx = 'AM_INV_TOTALT_Q_N';

%series_gnp = 'GNP.Q_N';
series_gnp = 'GNP.Q_FP';

% residential and nonresidential equipment and structures, capital stock
series_ns = {'AM_NSTOCK_TOTALT_Q_N'}; % includes inventories... not exactly like Fama
series_roc = 'ROC_BS.Q';
series_roc_denom = {'AM_NSTOCK_TOTALT_Q_N'}; % series_ns_fred + inventories
%{
**************************************
NET CAPITAL STOCK

from:
https://www.bea.gov/resources/learning-center/definitions-and-introduction-fixed-assets

With the perpetual inventory method that is used to derive the estimates 
presented here, the net stock in the historical-cost valuation and (at 
the individual asset) in the real-cost valuation, which are described 
below, is calculated as the cumulative value of past investment less the 
cumulative value of past depreciation.

Net stock in current-cost valuation is the value of the items in the
real-cost net stock measured in the prices of the current
end of year.7 

Average age of net stock for a given end of year is a weighted average of 
the ages of all investment in the stock at that yearend, with the weight 
for each age based on its value in the net stock.

7. The difference between gross and net stocks of an asset in current-cost 
valuation is usually not equal to the sum of past depreciation charges on
that asset in current-cost valuation because the difference between the 
two stock measures is a function of prices in only the current year while 
the past depreciation charges are a function of prices in all years since 
the asset was purchased.

For example, Residential structures (RCVSRNWMVBSNNCB, Source ID: FL105012665.Q):
         
https://www.federalreserve.gov/apps/fof/SeriesAnalyzer.aspx?s=FL105012665&t=
Last edited on: 09/17/2013

Series analyzer for FA105012665.Q
Year-end level from unpublished BEA data. Other quarter's levels are 
derived from year-end level, plus the quarterly unadjusted transactions 
calculated as gross investment ( FOF series FU105012063) less depreciation 
(FOF series FU106320063), plus an estimate of capital gains based on NIPA,
Table 5.3.4. Price Indexes for Private Fixed Investment by Type, 
line 21, Residential Structures.

**************************************
TAXES
Tax liability can be computed in one of two ways:
(1) subtracting profit after tax from profit before tax
(2) gauging it from total tax liability (BOGZ1FL143178005Q)

These are close but not exactly the same.



**************************************
RETURN ON CAPITAL, ROC
cash flow / net capital stock

Fama (1981) writes:
The numerator is the before-tax cash flow to nonfinancial corporations less
the BEA's estimate of replacement cost depreciation and less federal,
state, and local corporate profits tax liabilities.

cash flow =
before-tax cash flow to nonfinancial corporations
-
BEA's estimate of replacement cost depreciation
-
federal, state, and local corporate profits tax liabilities.

before-tax cash flow to nonfinancial corporations =
(1) profits before taxes 
(PBT, from IRS, includes depreciation for tax purposes in the negative)
+
(2) depreciation taken for tax purposes
+
(3) monetary net interest paid


My understanding is that CCAdj reflects 
CCAdj
=
BEA's estimate of replacement cost depreciation
-
depreciation taken for tax purposes.

I calculated CCAdj in the following way:
CCAdj = BOGZ1FU106060035Q - BOGZ1FU106060005Q - NCBIVDQ027S
PBT incl. CCAdj&IVA - PBT excl. CCAdj&IVA - IVA

This is different from B470RC1Q027SBEA. I think this is because my CCAdj
is based on transactions, whereas B470RC1Q027SBEA is somehow a level


My final cashflow is: (doesn't include IVA but does include CCAdj)
PBT + CCAdj - tax + interest piad

BOGZ1FU106060005Q + CCAdj - BOGZ1FL143178005Q + BOGZ1FA106130001Q




**************************************
CORPORATE PROFITS

** PBT—sometimes referred to as “book profits”—reflects corporate income
regardless of any redistribution of income made through taxes. The PBT
estimates are primarily based on tax-return information provided by the IRS
in Statistics of Income: Corporation Income Tax Returns, which is then
adjusted to conform to BEA coverage and definitions. PBT is distributed to
government as taxes on corporate income and to shareholders in other
sectors as dividends or is retained as undistributed profits. The estimates
of PBT are prior to the IVA and CCAdj, so PBT reflects the charges used in
**business tax accounting for inventory withdrawals and for depreciation.**

Taxes on corporate income consists of taxes paid on corporate earnings to
federal, state, and local governments and to foreign governments. These
earnings include capital gains and other income excluded from PBT. The
taxes are measured on an accrual basis, net of applicable tax credits.

Corporate profits with IVA is defined in the same way as corporate profits
with IVA and CCAdj, except corporate profits with IVA reflects the
depreciation-accounting practices used for federal income tax returns.
Profits by industry is shown on this basis because estimates of the CCAdj
by industry are not available.

Profits after tax without IVA and CCAdj is equal to PBT less taxes on
corporate income. It consists of net dividends and undistributed corporate
profits. This measure is often used in comparisons with the S&P measures of
reported earnings.

CCAdj. Depreciation is an expense in calculating profits. Depreciation
measured on a business-accounting basis must be adjusted to reflect
consistent economic-accounting measures that are valued at
current-replacement cost. The CCAdj is a two-part adjustment that:

(1) converts valuations of depreciation that are based on a mixture of
service lives and depreciation patterns specified in the tax code to
valuations that are based on uniform service lives and empirically based
depreciation patterns;
(2) like the IVA, converts the measures of depreciation to a current-cost
basis by removing from profits the gain or loss that arises from valuing
the depreciation of fixed assets at the prices of earlier periods.
%}

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
    },...
    'VariableNames', cell(1,14)...
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

    min_date = sample_tab.min_date(sample_tab.period==period);
    max_date = sample_tab.max_date(sample_tab.period==period);
    
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

    EIPF_ = T.EIPF/100*pcent_mult;
    UIPF_ = T.UIPF/100*pcent_mult;
    EIAR_ = T.EIAR1_101/100*pcent_mult;
    UIAR_ = T.UIAR1_101/100*pcent_mult;
    EITB_ = T.EITB/100*pcent_mult;
    UITB_ = T.UITB/100*pcent_mult;
    T = removevars(T, {'EITB', 'UITB', 'EIPF', 'UIPF', 'EIAR1_101', 'UIAR1_101'});  

    eval(sprintf('EI = EI%s_;', char(EI_measure)))
    eval(sprintf('UI = UI%s_;', char(EI_measure)))



    
    % ***************************************
    % Capital expenditures

    idx = ismember({data.series}, series_cx);
    CX = T{:, data(idx).VariableNames};

    % ***************************************
    % Net capital stock

    idx = ismember({data.series}, series_ns);
    NS = sum(T{:, [data(idx).VariableNames]}, 2);

    % ***************************************
    % Return on Capital, ROC
        
    idx = ismember({data.series}, series_roc);
    ROC = sum(T{:, [data(idx).VariableNames]}, 2);
        
    CX_NS = CX ./ NS;

    dROC = [NaN;...
        ROC(2:end)./...
        ROC(1:end-1)];    
    
    dCX_NS = [NaN;...
        CX_NS(2:end)./...
        CX_NS(1:end-1)];    
    
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

    idx = ismember({data.series}, series_CPI);
    CPI = T{:, data(idx).VariableNames};
    Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
   
    if log_approx

        Pi = log(Pi)*pcent_mult;
        dPR = log(dPR)*pcent_mult;
        dRGNP = log(dRGNP)*pcent_mult;
        dM = log(dM)*pcent_mult;
        dCX_NS = log(dCX_NS)*pcent_mult;
        dROC = log(dROC)*100;

    else

        Pi = (Pi - 1)*pcent_mult;
        dPR = (dPR - 1)*pcent_mult;
        dRGNP = (dRGNP- 1)*pcent_mult;
        dM = (dM - 1)*pcent_mult;
        dCX_NS = (dCX_NS - 1)*pcent_mult;
        dROC = (dROC- 1)*pcent_mult;

    end

    % define leads and lags
    % T = sortrows(T, 'date', 'ascending');  % make sure T is ascending
    dPR_lag = [NaN; dPR(1:end-1)];
    dM_lag = [NaN; dM(1:end-1)];
    dCX_NS_lag = [NaN; dCX_NS(1:end-1)];
    dROC_lag = [NaN; dROC(1:end-1)];
    dRGNP_lag = [NaN; dRGNP(1:end-1)];
    Rf_lag = [NaN; Rf(1:end-1)];

    dPR_lead = [dPR(2:end); NaN];
    dM_lead = [dM(2:end); NaN];
    dCX_NS_lead = [dCX_NS(2:end); NaN];
    dROC_lead = [dROC(2:end); NaN];
    dRGNP_lead = [dRGNP(2:end); NaN];
    Rf_lead = [Rf(2:end); NaN];
    
    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    EI = EI(idx_sample);
    UI = UI(idx_sample);

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
    
    dCX_NS = dCX_NS(idx_sample);
    dCX_NS_lead = dCX_NS_lead(idx_sample);
    dCX_NS_lag = dCX_NS_lag(idx_sample);

    dROC = dROC(idx_sample);
    dROC_lead = dROC_lead(idx_sample);
    dROC_lag = dROC_lag(idx_sample);
    
    N = size(Pi, 1);
    table2_a = ...
    [...
        corr(Pi, [dRGNP_lag, dRGNP, dRGNP_lead], 'rows', 'complete');...
        corr(Pi, [dPR_lag, dPR, dPR_lead], 'rows', 'complete');...
        corr(Pi, [dCX_NS_lag, dCX_NS, dCX_NS_lead], 'rows', 'complete');...
        corr(Pi, [dROC_lag, dROC, dROC_lead], 'rows', 'complete')...
    ];
    table2_b = ...
    [...
        corr(EI, [dRGNP_lag, dRGNP, dRGNP_lead], 'rows', 'complete');...
        corr(EI, [dPR_lag, dPR, dPR_lead], 'rows', 'complete');...
        corr(EI, [dCX_NS_lag, dCX_NS, dCX_NS_lead], 'rows', 'complete');...
        corr(EI, [dROC_lag, dROC, dROC_lead], 'rows', 'complete')...
    ];
    table2_c = ...
    [...
        corr(dM, [dRGNP_lag, dRGNP, dRGNP_lead], 'rows', 'complete');...
        corr(dM, [dPR_lag, dPR, dPR_lead], 'rows', 'complete');...
        corr(dM, [dCX_NS_lag, dCX_NS, dCX_NS_lead], 'rows', 'complete');...
        corr(dM, [dROC_lag, dROC, dROC_lead], 'rows', 'complete')...
    ];
    
    min_dates = repmat(min_date, 4, 1);
    max_dates = repmat(max_date, 4, 1);
    Ynames = {'DRGNP', 'DPR', 'DCX/NS', 'DROC'}';
    VariableNames = {'tau=-1', 'tau=0', 'tau=+1'};
    Xnames = {'I', 'I', 'I', 'I'}';
    table2_a = [...
        table(min_dates, max_dates, 'VariableNames', {'min_date', 'max_date'}),...
        cell2table(Xnames, 'VariableNames', {'X_t'}),...
        cell2table(Ynames, 'VariableNames', {'Y_t+tau'}),...
        array2table(table2_a, 'VariableNames', VariableNames)...
        ];
    
    %Xnames = {'EITB', 'EITB', 'EITB', 'EITB'}';
    Xnames = repmat(cellstr(['EI', char(EI_measure)]), 4, 1);
    table2_b = [...
        table(min_dates, max_dates, 'VariableNames', {'min_date', 'max_date'}),...
        cell2table(Xnames, 'VariableNames', {'X_t'}),...
        cell2table(Ynames, 'VariableNames', {'Y_t+tau'}),...
        array2table(table2_b, 'VariableNames', VariableNames)...
        ];
    
    Xnames = {'BG', 'BG', 'BG', 'BG'}';
    table2_c = [...
        table(min_dates, max_dates, 'VariableNames', {'min_date', 'max_date'}),...
        cell2table(Xnames, 'VariableNames', {'X_t'}),...
        cell2table(Ynames, 'VariableNames', {'Y_t+tau'}),...
        array2table(table2_c, 'VariableNames', VariableNames)...
        ];
    

    table2_ = [table2_a; table2_b; table2_c];
    table2 = [table2; table2_];
    
end



fpath_table2 = fullfile(root, 'Fama1981_il', sprintf('table2.%s.csv', freq));

if exist(fpath_table2, 'file') == 2 && false
    table2_ = readtable(fpath_table2);
    table2_.Properties.VariableNames = table2.Properties.VariableNames;
    table2 = [table2; table2_];
end
%}
%{
% save each iteration alone

fpath_table2 = fullfile(root, 'Fama1981', 'table2.csv');
counter = 0;
while exist(fpath_table2, 'file') == 2
    fpath_tab1 = fullfile(root, 'Fama1981', sprintf('table2_%d.csv', counter));
    counter = counter + 1;
end
        
%}
%tab1 = sortrows(tab1, {'min_date', 'max_date'});
writetable(table2, fpath_table2)
table2



%{

This script runs replicates 7 regressions in table 3, Fama (1981).

It does so with updated data, which is a bit different to the original.

For example, some series have had revisions. I did the best I could to look
at the original data, even going to the Survey of Current Business real
time publications to see some of the series Fama used. In the end I decided
against using those series, because results are not very different when
using the updated series.

TODO:
Some work is needed here to adjust frequencies. Some data are monthly, some
are quarterly, in the end I need Annual frequency.

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
lib_gs = fullfile(lib_capital_stock, 'gstock');  % gross capital stock

%lib_ea = fullfile(lib_israel, 'indprod');  % economic activity
lib_ea = fullfile(lib_israel);  % economic activity, monthly series

lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M', 'with_metadata');
%lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M');
lib_pop = fullfile(lib_israel, 'POP');
lib_consumption = fullfile(lib_israel, 'C');
%lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');

lib_stock_markets = fullfile(lib_data, 'stock_markets');
lib_SM = fullfile(lib_israel, 'SM', freq', 'with_metadata');

lib_EITB = fullfile(root, 'EITB_il');
lib_UITB = fullfile(root, 'UITB_il');
lib_EIMD = fullfile(root, 'EIMD_il');
lib_UIMD = fullfile(root, 'UIMD_il');



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
sample_tab = table(...
    [1:3]',...
    [datetime(2000,1,31), datetime(2000,1,31), datetime(2008,10,31),...
     ]',...
    [datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}


%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

table6_1 = [];
table6_2 = [];
table6_3 = [];
table6_4 = [];
table6_5 = [];
table6_6 = [];
table6_7 = [];


series_CPI = 'CP_SA.M';
series_CPI_NSA = 'CP_NSA.M';EI_measure = {'AR'};

if strcmp(freq, 'M')
    series_Rf = 'TSB_BAGR_MAKAM_01M.M'; % 1 month MAKAM, annualized
elseif strcmp(freq, 'Q')
    series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized
elseif strcmp(freq, 'A')
    series_Rf = 'TSB_BAGR_MAKAM_12M.M'; % 12 month MAKAM, annualized
else
    error('freq must be M, Q or A')
end

series_ms = 'ZA108.M'; % waiting for Elad from MOS, not incl coin
%series_ms = {'A346.M_E', 'ZA215.M'}; % ask Elad about RESNALNS, 

series_ea = 'mdd13.level.m';
%series_ea = 'CLS11.TPR_C.Q_SA_CHAINED'; % excl diamonds


% note: there are discrepancies between NIPA and FRED. use NIPA
%series_cx = 'AM_INV_TOTALT_Q_N';
series_cx = 'AM_INV_ISKI_Q_N';

%series_gnp = 'GNP.Q_N';
series_gnp = 'GNP.Q_FP';

% residential and nonresidential equipment and structures, capital stock
series_ns = {'AM_NSTOCK_ISKI_Q_N'}; % includes inventories... not exactly like Fama
series_gs = {'AM_NSTOCK_ISKI_Q_N'}; % includes inventories... not exactly like Fama
series_roc = 'ROC_BS.Q';
series_roc_denom = {'AM_NSTOCK_TOTALT_Q_N'}; % series_ns_fred + inventories

%series_tase = 'SPASTT01ILM661N';
series_tase = 'MD_M2.M';


series_C = {'C_@DUR.Q_FP_SA'}; % excl durable, should contain services
% not needed
%series_PD = {'PCEPINDG', 'PCEPIS'}; % I have the real series, don't really
%need this. (implicit deflator)
%series_pop = 'DEM_POP_AVE_CHN.Q';



%{
**************************************
NET CAPITAL STOCK

There are big differences between BEA contemporary series of the net
capital stock and the series published in April 1976 at the Survey of
Current Business (SCB): Fixed Nonresidential Business and Residential Capital in
the United States, 1925-75, that is cited by Fama, 1981. (and updates of it
in the years between 1976-1981.

For example, the first observation in the series, Dec-1953:
Table 2.—Current-Dollar Net Stocks of Fixed Nonresidential Business Capital,
by Major Industry Group and Legal Form of Organization, 1925-75
Corporate, Nonfinancial:
Equipment+Structures: 177.0 (84.9+92.1)

Current NIPA, Table 4.1. Current-Cost Net Stock of Private Nonresidential
Fixed Assets by Industry Group and Legal Form of Organization
Equipment+Structures(+intelectual property): 288.4 (86.7+189.2)
the value of structures is almost double.

series_ns_fred agree with current NIPA tables of fixed assets.

%}

%{
**************************************
CAPITAL EXPENDITURE

Moreover, in the SCB I can't find capital expenditure. This is not merely
the difference between one year to another, because of depreciation and
price changes.

The stock estimates are derived by the perpetual inventory method, which 
starts with investment flows and calculates gross capital stock for any 
given yearend by cumulating past investment flows and deducting discards.


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
%}


%{
**************************************
TAXES
Tax liability can be computed in one of two ways:
(1) subtracting profit after tax from profit before tax
(2) gauging it from total tax liability (A054RC1A027NBEA)

These are close but not exactly the same.


%}


%{
**************************************
CONSUMPTION OF FIXED CAPITAL
Tax liability can be computed in one of two ways:
(1) subtracting profit after tax from profit before tax
(2) gauging it from total tax liability (A054RC1A027NBEA)

These are close but not exactly the same. I'm taking the official series:
Nonfinancial Corporate Business; Capital Consumption Allowance, Transactions


%}

%{
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


%}

%{

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
    'AM_NSTOCK_ISKI_Q_N',...
    'AM_GSTOCK_ISKIM_Q_N',... % Business Sector; Equipment
    'AM_GSTOCK_ISKIS_Q_N',... % Business Sector; Structures
    'AM_GSTOCK_BILDRES_Q_N',...  % Construction Sector; Residential All
    'AM_GSTOCK_TOTALT_Q_N',...
    'AM_GSTOCK_ISKI_Q_N',...
    series_roc,...
    'SPASTT01ILM661N',...
    'MD_M2.M',...
    'MD_M137.M',...
    'MD_M142.M',...
    'MD_M143.M',...
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
    lib_ns,...
    lib_gs,...
    lib_gs,...
    lib_gs,...
    lib_gs,...
    lib_gs,...
    lib_roc,...
    lib_stock_markets,...
    lib_SM,...
    lib_SM,...
    lib_SM,...
    lib_SM,...
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
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'FAME',...
    'FRED',...
    'FAME',...
    'FAME',...
    'FAME',...
    'FAME',...
    },...
    'VariableNames', cell(1,25)...
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
        %{
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true, 'Type', 'left');
        %}
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true');
        %metadata = [metadata; metadata_];
    end
    

end

lags_ar = [1:4, 8, 12];
% In-sample UIAR
fpath_factor = fullfile(lib_matlab_TASE, 'arima_101_il.csv');
arima_101 = readtable(fpath_factor);
fpath_factor = fullfile(lib_matlab_TASE, 'arima_011_il.csv');
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
    % to end of quarter
    EITB.date = EITB.date + calmonths(3);
    EITB.date = EITB.date -1;

    EITB1 = table2timetable(EITB1);
    EITB1 = retime(EITB1, 'quarterly', 'mean');
    EITB1 = timetable2table(EITB1);
    % to end of quarter
    EITB1.date = EITB1.date + calmonths(3);
    EITB1.date = EITB1.date - 1;

    EIAR = table2timetable(EIAR);
    EIAR = retime(EIAR, 'quarterly', 'mean');
    EIAR = timetable2table(EIAR);
    % to end of quarter
    EIAR.date = EIAR.date + calmonths(3);
    EIAR.date = EIAR.date - 1;

    EIPF = table2timetable(EIPF);
    EIPF = retime(EIPF, 'quarterly', 'mean');
    EIPF = timetable2table(EIPF);
    % to end of quarter
    EIPF.date = EIPF.date + calmonths(3);
    EIPF.date = EIPF.date - 1;

    % Unexpected inflation
    UITB = table2timetable(UITB);
    UITB = retime(UITB, 'quarterly', 'mean');
    UITB = timetable2table(UITB);
    % to end of quarter
    UITB.date = UITB.date + calmonths(3);
    UITB.date = UITB.date - 1;

    UITB1 = table2timetable(UITB1);
    UITB1 = retime(UITB1, 'quarterly', 'mean');
    UITB1 = timetable2table(UITB1);
    % to end of quarter
    UITB1.date = UITB1.date + calmonths(3);
    UITB1.date = UITB1.date - 1;

    UIAR = table2timetable(UIAR);
    UIAR = retime(UIAR, 'quarterly', 'mean');
    UIAR = timetable2table(UIAR);
    % to end of quarter
    UIAR.date = UIAR.date + calmonths(3);
    UIAR.date = UIAR.date - 1;

    UIPF = table2timetable(UIPF);
    UIPF = retime(UIPF, 'quarterly', 'mean');
    UIPF = timetable2table(UIPF);
    % to end of quarter
    UIPF.date = UIPF.date + calmonths(3);
    UIPF.date = UIPF.date - 1;




elseif strcmp(freq, 'A')
    lags_ar = 1:4;
    T = table2timetable(T);
    T = retime(T, 'yearly', 'lastvalue');
    T = timetable2table(T);

    EITB = table2timetable(EITB);
    EITB = retime(EITB, 'yearly', 'mean');
    EITB = timetable2table(EITB);
    EITB.date = EITB.date + calmonths(12);
    EITB.date = EITB.date - 1;

    EITB1 = table2timetable(EITB1);
    EITB1 = retime(EITB1, 'yearly', 'mean');
    EITB1 = timetable2table(EITB1);
    EITB1.date = EITB1.date - calmonths(12);
    EITB1.date = EITB1.date - 1;

    EIAR = table2timetable(EIAR);
    EIAR = retime(EIAR, 'yearly', 'mean');
    EIAR = timetable2table(EIAR);
    EIAR.date = EIAR.date - calmonths(12);
    EIAR.date = EIAR.date - 1;

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
    UITB.date = UITB.date - calmonths(12);

    UITB1 = table2timetable(UITB1);
    UITB1 = retime(UITB1, 'yearly', 'mean');
    UITB1 = timetable2table(UITB1);
    UITB1.date = UITB1.date - 1;
    UITB1.date = UITB1.date - calmonths(12);

    UIAR = table2timetable(UIAR);
    UIAR = retime(UIAR, 'yearly', 'mean');
    UIAR = timetable2table(UIAR);
    UIAR.date = UIAR.date - 1;
    UIAR.date = UIAR.date - calmonths(12);

    UIPF = table2timetable(UIPF);
    UIPF = retime(UIPF, 'yearly', 'mean');
    UIPF = timetable2table(UIPF);
    UIPF.date = UIPF.date - 1;
    UIPF.date = UIPF.date - calmonths(12);



else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

T = fillmissing(T, 'previous');





% Sample period
%for period = sample_tab.period(end-2)'
for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);

    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;


    T = outerjoin(T, EIPF, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, UIPF, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, EIAR, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, UIAR, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);

    EIPF_ = T.EIPF/100;
    UIPF_ = T.UIPF/100;
    EIAR_ = T.EIAR1_101/100;
    UIAR_ = T.UIAR1_101/100;

    T = removevars(T, {'EIPF', 'UIPF'});  
    T = removevars(T, {'EIAR1_101', 'UIAR1_101'});  

    
    if strcmp(freq, 'A')
        EIAR_ = EIAR_*12;
        UIAR_ = UIAR_*12;
        EIPF_ = EIPF_*12;
        UIPF_ = UIPF_*12;
    elseif strcmp(freq, 'Q')
        EIAR_ = EIAR_*4;
        UIAR_ = UIAR_*4;
        EIPF_ = EIPF_*4;
        UIPF_ = UIPF_*4;
    
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
    
    CX_NS = CX ./ NS*pcent_mult;

    % ***************************************
    % Return on Capital, ROC
        
    idx = ismember({data.series}, series_roc);
    ROC = T{:, data(idx).VariableNames}*pcent_mult;
    
    % ***************************************
    % CRSP value weighted inedx

    idx = ismember({data.series}, series_tase);
    TASE = T{:, data(idx).VariableNames};
    
    % ***************************************
    % Risk free rate

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames}/100;
    T.Rf = T{:, data(idx).VariableNames}/100;

    if strcmp(freq, 'M')
        T.Rf = T.Rf/12;
    elseif strcmp(freq, 'Q')
        T.Rf = T.Rf/4;
    elseif strcmp(freq, 'A')
        % pass
    else
        error('freq must be M, Q or A')
    end

    % ***************************************
    % Real GNP
    
    idx = ismember({data.series}, series_gnp);
    RGNP = T{:, data(idx).VariableNames};
    dRGNP = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}]; 
    dRGNP3 = [NaN(3, 1);...
        T{4:end,   data(idx).VariableNames}./...
        T{1:end-3, data(idx).VariableNames}];
    dRGNP4 = [NaN(4, 1);...
        T{5:end,   data(idx).VariableNames}./...
        T{1:end-4, data(idx).VariableNames}];
    dRGNP12 = [NaN(12,1);...
        T{13:end,   data(idx).VariableNames}./...
        T{1:end-12, data(idx).VariableNames}];

    dRGNP_f = [dRGNP(2:end); NaN(1, 1)];
    dRGNP3_f = [dRGNP3(4:end); NaN(3, 1)];
    dRGNP4_f = [dRGNP4(5:end); NaN(4, 1)];
    dRGNP12_f = [dRGNP12(13:end); NaN(12, 1)];

    % ***************************************
    % Industrial production
    
    idx = ismember({data.series}, series_ea);
    PR = T{:, data(idx).VariableNames};
    dPR = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];

    dPR3 = [NaN(3, 1);...
        T{4:end,   data(idx).VariableNames}./...
        T{1:end-3, data(idx).VariableNames}];
    dPR4 = [NaN(4, 1);...
        T{5:end,   data(idx).VariableNames}./...
        T{1:end-4, data(idx).VariableNames}];
    dPR12 = [NaN(12,1);...
        T{13:end,   data(idx).VariableNames}./...
        T{1:end-12, data(idx).VariableNames}];

    dPR_f = [dPR(2:end); NaN(1, 1)];
    dPR3_f = [dPR3(4:end); NaN(3, 1)];
    dPR4_f = [dPR4(5:end); NaN(4, 1)];
    dPR12_f = [dPR12(13:end); NaN(12, 1)];


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

    EIAR_(isnan(EIAR_)) = Pi(isnan(EIAR_))-1;
    UIAR_(isnan(UIAR_)) = 0;
    
    
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
    dROC3 = [NaN(3, 1);...
        ROC(4:end)-...
        ROC(1:end-3)];
    dROC4 = [NaN(4, 1);...
        ROC(5:end)-...
        ROC(1:end-4)];
    dROC12 = [NaN(12,1);...
        ROC(13:end)-...
        ROC(1:end-12)];    
    dROC = [NaN;...
        ROC(2:end)-...
        ROC(1:end-1)];  

    dROC_f = [dROC(2:end); NaN(1, 1)];
    dROC3_f = [dROC3(4:end); NaN(3, 1)];
    dROC4_f = [dROC4(5:end); NaN(4, 1)];
    dROC12_f = [dROC12(13:end); NaN(12, 1)];

    dCX_NS = [NaN;...
        CX_NS(2:end)-...
        CX_NS(1:end-1)];  
    dCX_NS3 = [NaN(3, 1);...
        CX_NS(4:end)-...
        CX_NS(1:end-3)];
    dCX_NS4 = [NaN(4, 1);...
        CX_NS(5:end)-...
        CX_NS(1:end-4)];
    dCX_NS12 = [NaN(12,1);...
        CX_NS(13:end)-...
        CX_NS(1:end-12)];

    dCX_NS_f = [dCX_NS(2:end); NaN(1, 1)];
    dCX_NS3_f = [dCX_NS3(4:end); NaN(3, 1)];
    dCX_NS4_f = [dCX_NS4(5:end); NaN(4, 1)];
    dCX_NS12_f = [dCX_NS12(13:end); NaN(12, 1)];


    dTASE = [NaN;...
        TASE(2:end)./...
        TASE(1:end-1)];
    
    if log_approx
        
        Pi = log(Pi)*pcent_mult;

        dPR = log(dPR)*pcent_mult;
        dPR3 = log(dPR3)*pcent_mult;
        dPR4 = log(dPR4)*pcent_mult;
        dPR12 = log(dPR12)*pcent_mult;
        dPR3_f = log(dPR3_f)*pcent_mult;
        dPR12_f = log(dPR12_f)*pcent_mult;

        dRGNP = log(dRGNP)*pcent_mult;
        dRGNP3 = log(dRGNP3)*pcent_mult;
        dRGNP4 = log(dRGNP4)*pcent_mult;
        dRGNP12 = log(dRGNP12)*pcent_mult;
        dRGNP3_f = log(dRGNP3_f)*pcent_mult;
        dRGNP12_f = log(dRGNP12_f)*pcent_mult;

        dM = log(dM)*pcent_mult;
        dTASE = log(dTASE)*pcent_mult;
        %dCX_NS = log(dCX_NS)*pcent_mult;
        %dROC = log(dROC)*pcent_mult;

    else

        Pi = (Pi - 1)*pcent_mult;

        dPR = (dPR - 1)*pcent_mult;
        dPR3 = (dPR3 - 1)*pcent_mult;
        dPR4 = (dPR4 - 1)*pcent_mult;
        dPR12 = (dPR12 - 1)*pcent_mult;
        dPR3_f = (dPR3_f - 1)*pcent_mult;
        dPR12_f = (dPR12_f - 1)*pcent_mult;


        dRGNP = (dRGNP- 1)*pcent_mult;
        dRGNP3 = (dRGNP- 1)*pcent_mult;
        dRGNP4 = (dRGNP- 1)*pcent_mult;
        dRGNP12 = (dRGNP- 1)*pcent_mult;
        dRGNP3_f = (dRGNP- 1)*pcent_mult;
        dRGNP12_f = (dRGNP- 1)*pcent_mult;

        dM = (dM - 1)*pcent_mult;
        dTASE = (dTASE- 1)*pcent_mult;
        %dCX_NS = (dCX_NS - 1)*pcent_mult;
        %dROC = (dROC- 1)*pcent_mult;

    end
    
    RS = dTASE-Pi; % real stock return
    XR = dTASE-Rf;
%{
From Fama, 1981, footnote 4.
The real stock return RS, is the annual continuously compounded nominal 
return on a value weighted portfolio of all New York Stock Exchange common
stocks less the annual continuously comnpounded inflation rate calculated
from the U.S. Consumer Price Index. The nominal stock return is from the 
Center for Research in Security Prices of the University of Chicago. The 
monthly and quarterly versions of RS, used in later sections are defined 
similarly.
    
%}
    
    EIPF_ = EIPF_*pcent_mult;
    UIPF_ = UIPF_*pcent_mult;
    EIAR_ = EIAR_*pcent_mult;
    UIAR_ = UIAR_*pcent_mult;

 
    % define leads and lags
    % T = sortrows(T, 'date', 'ascending');  % make sure T is ascending
    dPR_lag = [NaN; dPR(1:end-1)];
    dM_lag = [NaN; dM(1:end-1)];
    dCX_NS_lag = [NaN; dCX_NS(1:end-1)];
    dROC_lag = [NaN; dROC(1:end-1)];
    dRGNP_lag = [NaN; dRGNP(1:end-1)];
    Rf_lag = [NaN; Rf(1:end-1)];
    EIPF_lag = [NaN; EIPF_(1:end-1)];
    UIPF_lag = [NaN; UIPF_(1:end-1)];
    EIAR_lag = [NaN; EIAR_(1:end-1)];
    UIAR_lag = [NaN; UIAR_(1:end-1)];

    dPR_lead = [dPR(2:end); NaN];
    dM_lead = [dM(2:end); NaN];
    dCX_NS_lead = [dCX_NS(2:end); NaN];
    dROC_lead = [dROC(2:end); NaN];
    dRGNP_lead = [dRGNP(2:end); NaN];
    Rf_lead = [Rf(2:end); NaN];

    CPI = CPI(idx_sample);
    Pi = Pi(idx_sample);
    Rf = Rf(idx_sample);
    
    dM = dM(idx_sample);
    
    dPR = dPR(idx_sample);
    dPR3 = dPR3(idx_sample);
    dPR4 = dPR4(idx_sample);
    dPR12 = dPR12(idx_sample);
    dPR_f = dPR_f(idx_sample);
    dPR3_f = dPR3_f(idx_sample);
    dPR4_f = dPR4_f(idx_sample);
    dPR12_f = dPR12_f(idx_sample);
    dPR_lead = dPR_lead(idx_sample);
    dPR_lag = dPR_lag(idx_sample);

    dRGNP = dRGNP(idx_sample);
    dRGNP3 = dRGNP3(idx_sample);
    dRGNP4 = dRGNP4(idx_sample);
    dRGNP12 = dRGNP12(idx_sample);
    dRGNP_f = dRGNP_f(idx_sample);
    dRGNP3_f = dRGNP3_f(idx_sample);
    dRGNP4_f = dRGNP4_f(idx_sample);
    dRGNP12_f = dRGNP12_f(idx_sample);
    dRGNP_lead = dRGNP_lead(idx_sample);
    dRGNP_lag = dRGNP_lag(idx_sample);
    
    CX_NS = CX_NS(idx_sample);
    dCX_NS = dCX_NS(idx_sample);
    dCX_NS_lead = dCX_NS_lead(idx_sample);
    dCX_NS_lag = dCX_NS_lag(idx_sample);

    ROC = ROC(idx_sample);
    dROC = dROC(idx_sample);
    dROC_lead = dROC_lead(idx_sample);
    dROC_lag = dROC_lag(idx_sample);
    
    TASE = TASE(idx_sample);
    dTASE = dTASE(idx_sample);
    RS = RS(idx_sample);
    XR = XR(idx_sample);

    EIPF_ = EIPF_(idx_sample);
    EIPF_lag = EIPF_lag(idx_sample);
    UIPF_ = UIPF_(idx_sample);
    UIPF_lag = UIPF_lag(idx_sample);
    EIAR_ = EIAR_(idx_sample);
    EIAR_lag = EIAR_lag(idx_sample);
    UIAR_ = UIAR_(idx_sample);
    UIAR_lag = UIAR_lag(idx_sample);
    dates = T.date(idx_sample);
    
    N = size(Pi, 1);

    if strcmp(freq, 'M')
        dPR_f_X = dPR12_f;
        dRGNP_f_X = dRGNP12_f;
        dROC_f_X = dROC12_f;
        dCX_NS_f_X = dCX_NS12_f;
    elseif strcmp(freq, 'Q')
        dPR_f_X = dPR4_f;
        dRGNP_f_X = dRGNP4_f;
        dROC_f_X = dROC4_f;
        dCX_NS_f_X = dCX_NS4_f;
    elseif strcmp(freq, 'A')
        dPR_f_X = dPR_f;
        dRGNP_f_X = dRGNP_f;
        dROC_f_X = dROC_f;
        dCX_NS_f_X = dCX_NS_f;
    end

    %{
    % code to check expected inflation and inflation are denominated in the
    % same units
    clf
    subplot(1,2,1)
    hold on
    plot(dates, Pi)
    plot(dates, EIAR_)
    plot(dates, EIPF_)
    subplot(1,2,2)
    hold on
    plot(dates, RS)
    plot(dates, XR)
    return
    %}
    
    % **********************************
    % REGRESSION (1)
    Y = RS;
    X = [ones(N, 1), EIAR_, UIAR_];

    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    table6_1_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), lags_ar)];
    table6_1 = [table6_1; table6_1_];

    
    % **********************************
    % REGRESSION (2)
    Y = RS;
    X = [ones(N, 1), EIPF_, UIPF_];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    table6_2_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), lags_ar)];
    table6_2 = [table6_2; table6_2_];


    % **********************************
    % removed dM (Base growth of money)
    % replaced PR with RGNP
    % REGRESSION (3)
    Y = RS;
    X = [ones(N, 1), dRGNP, dRGNP_f_X, UIPF_];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    table6_3_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), lags_ar)];
    table6_3 = [table6_3; table6_3_];

  
    % **********************************
    % removed dM (Base growth of money)
    % replaced PR with RGNP
    % REGRESSION (4)
    Y = RS;
    X = [ones(N, 1), dRGNP_f_X, EIPF_, UIPF_];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    table6_4_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), lags_ar)];
    table6_4 = [table6_4; table6_4_];


    % **********************************
    % removed dM (Base growth of money)
    % replaced PR with RGNP
    % REGRESSION (5)
    Y = RS;
    X = [ones(N, 1), dRGNP_f_X, EIAR_, UIAR_];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    table6_5_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), lags_ar)];
    table6_5 = [table6_5; table6_5_];


    % **********************************
    % REGRESSION (6)
    Y = RS;
    X = [ones(N, 1), dPR_f_X, EIAR_, UIAR_];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    table6_6_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), lags_ar)];
    table6_6 = [table6_6; table6_6_];


    % **********************************
    % replaced PR with RGNP
    % REGRESSION (7)
    Y = RS;
    X = [ones(N, 1), dRGNP_f_X, dM, UIAR_];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    table6_7_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), lags_ar)];
    table6_7 = [table6_7; table6_7_];


    
end



rhos = strcat('rho_', arrayfun(@num2str, lags_ar, 'UniformOutput', false));



table6_1 = array2table(table6_1, 'VariableNames', [{...
    'nobs', 'constant', 'EIAR_{t-1}', 'UIAR_t',...
    'tstat_constant', 'tstat_EIAR_{t-1}', 'tstat_UIAR_t',...
    'R_sq', 'R_sq-adj', 's(epsilon)'}, rhos...
    ]);

table6_2 = array2table(table6_2, 'VariableNames', [{...
    'nobs', 'constant', 'EIPF_t', 'UIPF_t',...
    'tstat_constant', 'tstat_EIPF_t', 'tstat_UIPF_t',...
    'R_sq', 'R_sq-adj', 's(epsilon)'}, rhos...
    ]);

table6_3 = array2table(table6_3, 'VariableNames', [{...
    'nobs', 'constant', 'DRGNP_t', 'DRGNP_{t+12}', 'UIPF_t',...
    'tstat_constant', 'tstat_DRGNP_t', 'tstat_DRGNP_{t+12}',...
    'tstat_UIPF_t',...
    'R_sq', 'R_sq-adj', 's(epsilon)'}, rhos...
    ]);

table6_4 = array2table(table6_4, 'VariableNames', [{...
    'nobs', 'constant', 'DRGNP_{t+12}', 'EIPF_t', 'UIPF_t',...
    'tstat_constant', 'tstat_DRGNP_{t+12}',...
    'tstat_EIPF_t', 'tstat_UIPF_t',...
    'R_sq', 'R_sq-adj', 's(epsilon)'}, rhos...
    ]);


table6_5 = array2table(table6_5, 'VariableNames', [{...
    'nobs', 'constant', 'DRGNP_{t+12}', 'EIAR_{t-1}', 'UIAR_t',...
    'tstat_constant', 'tstat_DRGNP_{t+12}',...
    'tstat_EIAR_{t-1}', 'tstat_UIAR_t',...
    'R_sq', 'R_sq-adj', 's(epsilon)'}, rhos...
    ]);


table6_6 = array2table(table6_6, 'VariableNames', [{...
    'nobs', 'constant', 'DPR_{t+12}', 'EIAR_{t-1}', 'UIAR_t',...
    'tstat_constant', 'tstat_DPR_{t+12}',...
    'tstat_EIAR_{t-1}', 'tstat_UIAR_t',...
    'R_sq', 'R_sq-adj', 's(epsilon)'}, rhos...
    ]);

table6_7 = array2table(table6_7, 'VariableNames', [{...
    'nobs', 'constant', 'DRGNP_{t+12}', 'BG_t', 'UIAR_t',...
    'tstat_constant', 'tstat_DRGNP_{t+12}',...
    'tstat_BG_t', 'tstat_UIAR_t',...
    'R_sq', 'R_sq-adj', 's(epsilon)'}, rhos...
    ]);




for tab_num = 1:7
    
fpath_table6 = fullfile(root, 'Fama1981_il', sprintf('table6_%d.%s.csv', tab_num, freq));
eval(sprintf('tab = table6_%d;', tab_num))
tab = [sample_tab, tab];
writetable(tab, fpath_table6)
tab
    
end

%%
table6_4



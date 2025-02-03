%{

This script runs replicates desc stats in table 4, Fama (1981).

It does so with updated data, which is a bit different to the original,
and also with vintage data from the Survey of Current Business (SCB) of the
Bureau of Economic Analysis (BEA).

For example, some series have had revisions. I did the best I could to look
at the original data, even going to the Survey of Current Business real
time publications to see some of the series Fama used. In the end I decided
against using those series, because results are not very different when
using the updated series.

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

%lib_ir = fullfile(lib_israel);
lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M', 'with_metadata');
lib_pop = fullfile(lib_israel, 'POP');
lib_consumption = fullfile(lib_israel, 'C');
%lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');

lib_stock_markets = fullfile(lib_data, 'stock_markets');

lib_EITB = fullfile(root, 'EITB_il');
lib_UITB = fullfile(root, 'UITB_il');


% flag to override period in downstream scripts: false by default
lag_inflation = 0;
log_approx = true;
%log_approx = false;
freq = 'A';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
vintage = true;

pcent_mult = 1; % set as 100 for percent

sample_tab = table(...
    [1:3]',...
    [datetime(1999,12,31), datetime(1999,12,31), datetime(2008,10,31),...
     ]',...
    [datetime(2019,12,31), datetime(2008,9,30), datetime(2019,12,31),...
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

table4 = table;

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

series_tase = 'SPASTT01ILM661N';

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
    'AM_NSTOCK_ISKI_Q_N',...
    'AM_GSTOCK_ISKIM_Q_N',... % Business Sector; Equipment
    'AM_GSTOCK_ISKIS_Q_N',... % Business Sector; Structures
    'AM_GSTOCK_BILDRES_Q_N',...  % Construction Sector; Residential All
    'AM_GSTOCK_TOTALT_Q_N',...
    'AM_GSTOCK_ISKI_Q_N',...
    series_roc,...
    'SPASTT01ILM661N',...
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
    },...
    'VariableNames', cell(1,21)...
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
        %metadata = [metadata; metadata_(:, metadata.Properties.VariableNames)];
        
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
%for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);

    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;    

    
    %{
    % In-sample EITB
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
    NS = T{:, [data(idx).VariableNames]};    
    
    % ***************************************
    % Gross capital stock
    
    idx = ismember({data.series}, series_gs);
    GS = sum(T{:, [data(idx).VariableNames]}, 2); % current cost, gross, net; not summing, ipd...
    
    % first difference in gross stock = gross investment flows
    dGS = [NaN(1, size(GS,2));...
        GS(2:end, :)-...
        GS(1:end-1, :)];
    

    % ***************************************
    % Return on Capital, ROC
        
    idx = ismember({data.series}, series_roc);
    ROC = sum(T{:, [data(idx).VariableNames]}, 2);
    
    % ***************************************
    % CRSP value weighted inedx

    idx = ismember({data.series}, series_tase);
    TASE = T{:, data(idx).VariableNames};
    
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
        
    CX_NS = CX ./ sum(NS, 2);
        
    dROC = [NaN;...
        ROC(2:end, :)-...
        ROC(1:end-1, :)];
    
    
    dCX_NS = [NaN;...
        CX_NS(2:end, :)-...
        CX_NS(1:end-1, :)];    

    dCRSP = [NaN;...
        TASE(2:end, :)./...
        TASE(1:end-1, :)];
    

   
    if log_approx
        
        %Pi = log(Pi)*pcent_mult;
        dPR = log(dPR)*pcent_mult;
        dRGNP = log(dRGNP)*pcent_mult;
        dM = log(dM)*pcent_mult;
        dCRSP = log(dCRSP)*pcent_mult;
        %dCX_NS = log(dCX_NS)*pcent_mult;
        %dROC = log(dROC)*pcent_mult;

    else

        %Pi = (Pi - 1)*pcent_mult;
        dPR = (dPR - 1)*pcent_mult;
        dRGNP = (dRGNP- 1)*pcent_mult;
        dM = (dM - 1)*pcent_mult;
        dCRSP = (dCRSP- 1)*pcent_mult;
        %dCX_NS = (dCX_NS - 1)*pcent_mult;
        %dROC = (dROC- 1)*pcent_mult;

    end
    
    Pi = (Pi - 1)*pcent_mult;
    RS = dCRSP-Pi; % real stock return
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
    
    EI = EI/100*pcent_mult;
    UI = UI/100*pcent_mult;

    CPI = CPI(idx_sample);
    Pi = Pi(idx_sample);
    Rf = Rf(idx_sample);
    
    dM = dM(idx_sample);
    
    dPR = dPR(idx_sample);

    dRGNP = dRGNP(idx_sample);
    
    CX_NS = CX_NS(idx_sample);
    ROC = ROC(idx_sample);

    dCX_NS = dCX_NS(idx_sample);
    dROC = dROC(idx_sample);
    dCRSP = dCRSP(idx_sample);
    RS = RS(idx_sample);

    EI = EI(idx_sample);
    UI = UI(idx_sample);
    %{
    clf
    subplot(1,2,1)
    hold on
    %plot(T.date(idx_sample), CX_NS(idx_sample))
    idx = ismember({data.series}, series_ns);
    area(T.date(idx_sample), T{idx_sample, [data(idx).VariableNames]})
    plot(T.date(idx_sample), GS(idx_sample))
    T(idx_sample, [{'date'}, data(idx).VariableNames])
    yyaxis right
    plot(T.date(idx_sample), dGS(idx_sample))
    subplot(1,2,2)
    plot(T.date(idx_sample), CX_NS)
    
    sgtitle(sprintf('E[CX / NS]=%.5f; std[CX / NS]=%.5f',...
        mean(CX_NS, 'omitnan'), std(CX_NS, 'omitnan')))
    return
    %}
    
    N = size(Pi, 1);
    %{
    [mean(CX_NS, 'omitnan'), std(CX_NS, 'omitnan');...
        mean(dCX_NS, 'omitnan'), std(dCX_NS, 'omitnan')]
    return
    %}
    table4_ = [...
        arrayfun(@(lag) autocorr_(Pi, lag), [1:5]), mean(Pi, 'omitnan'), std(Pi, 'omitnan');...
        arrayfun(@(lag) autocorr_(dM, lag), [1:5]), mean(dM, 'omitnan'), std(dM, 'omitnan');...
        arrayfun(@(lag) autocorr_(CX_NS, lag), [1:5]), mean(CX_NS, 'omitnan'), std(CX_NS, 'omitnan');...
        arrayfun(@(lag) autocorr_(ROC, lag), [1:5]), mean(ROC, 'omitnan'), std(ROC, 'omitnan');...
        arrayfun(@(lag) autocorr_(dCX_NS, lag), [1:5]), mean(dCX_NS, 'omitnan'), std(dCX_NS, 'omitnan');...
        arrayfun(@(lag) autocorr_(dROC, lag), [1:5]), mean(dROC, 'omitnan'), std(dROC, 'omitnan');...
        arrayfun(@(lag) autocorr_(dRGNP, lag), [1:5]), mean(dRGNP, 'omitnan'), std(dRGNP, 'omitnan');...
        arrayfun(@(lag) autocorr_(dPR, lag), [1:5]), mean(dPR, 'omitnan'), std(dPR, 'omitnan');...
        arrayfun(@(lag) autocorr_(RS, lag), [1:5]), mean(RS, 'omitnan'), std(RS, 'omitnan');...
        ];
    
    VariableNames = {'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5', 'mean', 'std'};
    XNames = {'I_t', 'BG_t', 'CX_NS_t', 'ROC_t', 'dCX_NS_t', 'dROC_t', 'dRGNP_t', 'dPR_t', 'RS_t'}';
    table4_ = array2table(table4_, 'VariableNames', VariableNames);
    
    
    table4_ = [cell2table(XNames, 'VariableNames', {'Variable (x)'}), table4_];
    table4_ = [array2table(repmat(max_date, size(table4_, 1), 1), 'VariableNames', {'max_date'}), table4_];
    table4_ = [array2table(repmat(min_date, size(table4_, 1), 1), 'VariableNames', {'min_date'}), table4_];
    table4 = [table4; table4_];
    

    
end

table4
fpath_table4 = fullfile(root, 'Fama1981_il', sprintf('table4.%s.csv', freq));
writetable(table4, fpath_table4)




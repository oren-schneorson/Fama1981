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
lib_data = fullfile('/media', username, 'D', 'data');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')


lib_cpi = fullfile(lib_data, 'cpi_usa');
lib_ms = fullfile(lib_data, 'money_supply_usa');
lib_gnp = fullfile(lib_data, 'gnp');  % gross national product
lib_gnp_vintage = fullfile(lib_gnp, 'vintage_1978_07');  % gross national product
%lib_gnp_vintage = fullfile(lib_gnp, 'vintage_1980_01');  % gross national product
lib_ea = fullfile(lib_data, 'economic_activity');  % economic activity
lib_cx = fullfile(lib_data, 'cx');  % capital expenditure
lib_tax = fullfile(lib_data, 'tax');  % capital expenditure
lib_ns = fullfile(lib_data, 'ns', 'level');  % capital stock
%lib_ns_vintage = fullfile(lib_data, 'ns', 'vintage_1981_02');  % capital stock
lib_ns_vintage = fullfile(lib_data, 'ns', 'vintage_1979_08');  % capital stock
%lib_ns_vintage = fullfile(lib_data, 'ns', 'vintage_1978_09');  % capital stock
%lib_ns_vintage = fullfile(lib_data, 'ns', 'vintage_1977_08');  % capital stock
%lib_ns_vintage = fullfile(lib_data, 'ns', 'vintage_1976_08');  % capital stock
%lib_ns_vintage = fullfile(lib_data, 'ns', 'vintage_1976_04');  % capital stock

lib_ns_vintage_gross = lib_ns_vintage;


lib_cc = fullfile(lib_data, 'consumption_fixed_capital');  % capital consumption
lib_roc = fullfile(lib_data, 'roc');  % items for ROC_t
lib_mip = fullfile(lib_data, 'interest_paid');  % interest paid
lib_profits = fullfile(lib_data, 'profits');  % interest paid
lib_fi = fullfile(lib_data, 'fixed_investment');  % Fixed investment (gross)

lib_EITB = fullfile(root, 'EITB');
lib_UITB = fullfile(root, 'UITB');

lib_crsp = fullfile(lib_data, 'fama_french_factors', 'idx_mth');
lib_ir = fullfile(lib_data, 'ir');
lib_pop = fullfile(lib_data, 'population');
lib_consumption = fullfile(lib_data, 'consumption_usa');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');



% flag to override period in downstream scripts: false by default
lag_inflation = 0;
pcent_mult = 1; % set as 100 for percent
log_approx = true;
%log_approx = false;
freq = 'A';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
vintage = true;


if ~vintage
    lib_ns_vintage = lib_ns;
end


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
    [1:8]',...
    [datetime(1953,1,31), datetime(1953,1,31), datetime(2000,1,31),...
     datetime(1953,1,31), datetime(2007,7,31), datetime(1953,12,31),...
     datetime(1954,1,31), datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(2019,12,31), datetime(1977,12,31),...
    datetime(1976,12,31), datetime(1977,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});

%}

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

tab4 = table;


series_CPI = 'CPIAUCSL';
%series_CPI = 'CWUR0000SA0';

series_Rf = 'RF_FF';
series_crsp = 'CRSP_vw';
%series_ms = 'BOGMBASE';
series_ms = {'CURRCIR', 'RESBALNS'};
series_ea = 'INDPRO';

if strcmp(freq, 'A')
    series_pbt = 'A464RC1A027NBEA'; % profit before tax
    series_tax = 'BOGZ1FL103178005A'; % Nonfinancial Corporate Business; Total Taxes Payable; Liability, Level
    series_mip = 'B1101C1A027NBEA'; % Monetary interest paid: Domestic: Corporate business: Nonfinancial
    series_cc = 'BOGZ1FU106300015A'; % Nonfinancial Corporate Business; Capital Consumption Allowance, Transactions
    series_inv = 'NCBICBA027N'; % Nonfinancial Corporate Business; Inventories Excluding IVA, Current Cost Basis, Level
    series_cx_fred = 'BOGZ1FU105050005A'; % Nonfinancial Corporate Business; Total Capital Expenditures, Transactions
    %series_cx = {'BOGZ1FA105013063A', 'BOGZ1FA105013023A', 'BOGZ1FA105012063A', 'BOGZ1FA105012023A'}; % with breakdown of residential to equp & 
    series_cx = {'BOGZ1FA105013063A', 'BOGZ1FA105013023A', 'BOGZ1FA105012005A'}; % without breakdown of residential to equp & struct
    %series_gnp = 'A791RX0Q048SBEA'; % Real gross national product per capita
    series_gnp = 'GNPCA'; % Real Gross National Product

    % residential and nonresidential equipment and structures, capital stock
    %series_ns_fred = {'BOGZ1FL105013665A', 'BOGZ1FL105013265A',
    %'BOGZ1FL105012665A', 'BOGZ1FL105012265A'};  % annual series with breakdown to residential equipment and structures
    series_ns_fred = {'BOGZ1FL105013665A', 'BOGZ1FL105013265A', 'BOGZ1FL105012865A'}; % annual series, current cost (residential is equipment & struct)
    series_ns_fred_h = {'BOGZ1FL105013613A', 'BOGZ1FL105013213A', 'BOGZ1FL105012813A'}; % annual series, historical cost (residential is equipment & struct); vintage is constant dollar
    %series_ns_fred = {'BOGZ1FL105013665A', 'BOGZ1FL105013265A', 'BOGZ1FL105012665A'}; % annual series, current cost (residential is only struct)
    %series_ns_fred_h = {'BOGZ1FL105013613A', 'BOGZ1FL105013213A',
    %'BOGZ1FL105012613A'}; % annual series, historical cost (residential is only struct); no vintage
    series_roc_denom = [series_ns_fred, {'NCBICBA027N'}]; % series_ns_fred + inventories

elseif strcmp(freq, 'Q')
    series_pbt = 'A464RC1Q027SBEA'; % profit before tax
    series_tax = 'BOGZ1FL103178005Q'; % Nonfinancial Corporate Business; Total Taxes Payable; Liability, Level
    series_mip = 'BOGZ1FA106130001Q'; % Nonfinancial Corporate Business; Interest Paid, Transactions
    series_cc = 'BOGZ1FU106300015Q'; % Nonfinancial Corporate Business; Capital Consumption Allowance, Transactions
    series_inv = 'IABSNNCB'; % Nonfinancial Corporate Business; Inventories Excluding IVA, Current Cost Basis, Level
    series_cx_fred = 'BOGZ1FA105050005Q'; % Nonfinancial Corporate Business; Total Capital Expenditures, Transactions
    %series_cx = {'BOGZ1FA105013063Q', 'BOGZ1FA105013023Q', 'BOGZ1FA105012063Q', 'BOGZ1FA105012023Q'};
    series_cx = {'BOGZ1FA105013063Q', 'BOGZ1FA105013023Q', 'BOGZ1FA105012005Q'}; % residential without breakdown to equip & struct
    %series_gnp = 'A791RX0Q048SBEA'; % Real gross national product per capita
    series_gnp = 'GNPC96'; % Real Gross National Product
    
    % residential and nonresidential equipment and structures, capital stock
    %series_ns_fred = {'RCSNNWMVBSNNCB', 'BOGZ1FL105013265Q', 'RCVSRNWMVBSNNCB', 'BOGZ1FL105012265Q'}; % quarterly series with breakdown to residential equipment and structures
    series_ns_fred = {'RCSNNWMVBSNNCB', 'BOGZ1FL105013265Q', 'BOGZ1FL105012865Q'}; % quarterly series
    series_roc_denom = [series_ns_fred, {'IABSNNCB'}]; % series_ns_fred + inventories


else
    error('freq must be Q or A.') 
end


% note: there are discrepancies between NIPA and FRED. use NIPA
series_cx_nipa = 'BOGZ1FA105050005A_NIPA';

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

series_C = {'PCENDG', 'PCES'};
series_PD = {'PCEPINDG', 'PCEPIS'};


data = struct(...
    'series', {...
    'CWUR0000SA0',... CPI, not seasonally adjusted, BLS
    'CPIAUCSL',... Price index of all urban consumers
    series_crsp,...
    series_Rf,...
    'BOGMBASE',...
    'CURRCIR',...
    'RESBALNS',...
    'GNPCA',... % Annual real GNP
    'GNPC96',... % Quarterly real GNP
    'A791RX0Q048SBEA',... % real GNP per capita
    series_ea,...
    series_cx_nipa,...
    series_cx_fred,...
    'BOGZ1FA105013063A',... % Nonfinancial Corporate Business, Structures, nonresidential, annual freq
    'BOGZ1FA105013023A',... % Nonfinancial Corporate Business, Equipment, nonresidential, annual freq
    'BOGZ1FA105012063A',... % Nonfinancial Corporate Business, Structures, residential (incl. REIT), annual freq
    'BOGZ1FA105012023A',... % Nonfinancial Corporate Business, Equipment, residential, annual freq
    'BOGZ1FA105013063Q',... % Nonfinancial Corporate Business, Structures, nonresidential, Quarterly freq
    'BOGZ1FA105013023Q',... % Nonfinancial Corporate Business, Equipment, nonresidential, Quarterly freq
    'BOGZ1FA105012063Q',... % Nonfinancial Corporate Business, Structures, residential (incl. REIT), Quarterly freq
    'BOGZ1FA105012023Q',... % Nonfinancial Corporate Business, Equipment, residential, Quarterly freq
    'BOGZ1FL105013665A',... % Nonfinancial Corporate Business; Nonresidential Structures, Current Cost Basis, Level (vintage)
    'BOGZ1FL105013265A',... % Nonfinancial Corporate Business; Nonresidential Equipment, Current Cost Basis, Level  (vintage)
    'BOGZ1FL105012865A',... % Nonfinancial Corporate Business; Residential Equipment and Structures, Current Cost Basis, Level, annual freq  (vintage)
    'BOGZ1FL105013665A_gross',... % Nonfinancial Corporate Business; Nonresidential Structures, Current Cost Basis, Gross Level, annual freq (vintage)
    'BOGZ1FL105013265A_gross',... % Nonfinancial Corporate Business; Nonresidential Equipment, Current Cost Basis, Gross Level, annual freq (vintage)
    'BOGZ1FL105012865A_gross',... % Nonfinancial Corporate Business; Residential Equipment and Structures, Current Cost Basis, Gross Level, annual freq (vintage)
    'BOGZ1FL105013613A',... % Nonfinancial Corporate Business; Nonresidential Structures, Historical Cost Basis, Level (vintage)
    'BOGZ1FL105013213A',... % Nonfinancial Corporate Business; Nonresidential Equipment, Historical Cost Basis, Level (vintage)
    'BOGZ1FL105012813A',... % Nonfinancial Corporate Business; Residential Equipment and Structures, Historical Cost Basis, Level (vintage)
    'BOGZ1FL105013613A_gross',... % Nonfinancial Corporate Business; Nonresidential Structures, Historical Cost Basis, Gross Level (vintage)
    'BOGZ1FL105013213A_gross',... % Nonfinancial Corporate Business; Nonresidential Equipment, Historical Cost Basis, Gross Level (vintage)
    'BOGZ1FL105012813A_gross',... % Nonfinancial Corporate Business; Residential Equipment and Structures, Historical Cost Basis, Gross Level (vintage)
    'IABSNNCB',... % Nonfinancial Corporate Business; Inventories Excluding IVA
    'RCSNNWMVBSNNCB',... % Nonfinancial Corporate Business; Nonresidential Structures, Current Cost Basis, Level
    'BOGZ1FL105013265Q',... % Nonfinancial Corporate Business; Nonresidential Equipment, Current Cost Basis, Level
    'RCVSRNWMVBSNNCB',...  % Nonfinancial Corporate Business; Residential Structures, Current Cost Basis, Level
    'BOGZ1FL105012265Q',...  % Nonfinancial Corporate Business; Residential Equipment, Current Cost Basis, Level
    'BOGZ1FL105012665A',...  % Nonfinancial Corporate Business; Residential Structures, Current Cost Basis, Level
    'BOGZ1FL105012265A',...  % Nonfinancial Corporate Business; Residential Equipment, Current Cost Basis, Level
    'BOGZ1FL105012865Q',... % Nonfinancial Corporate Business; Residential Equipment and Structures, Current Cost Basis, Level, quarterly freq
    'BOGZ1FL105012213A',... % Nonfinancial Corporate Business; Residential Equipment, Historical Cost Basis, Level
    'BOGZ1FL105012613A',... % Nonfinancial Corporate Business; Residential Structures, Historical Cost Basis, Level
    'A054RC1A027NBEA',... % Taxes on corporate income, NIPAs, Billions of Dollars, Annual, Not Seasonally Adjusted
    series_tax,... % Nonfinancial Corporate Business; Total Taxes Payable; Liability, Level
    'BOGZ1FL143178005Q',... % Nonfinancial corporate business: Total tax liability
    'A053RC1A027NBEA',... % Corporate profits: Profits before taxes, NIPAs
    series_pbt,... % Nonfinancial corporate business: Profits before tax (without IVA and CCAdj), ְֳּ
    'A677RC1A027NBEA',... % Corporations: Capital consumption allowances, NIPAs, 
    series_cc... % Nonfinancial Corporate Business; Capital Consumption Allowance, Transactions
    series_mip,... % Nonfinancial Corporate Business; Interest Paid
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_crsp,...
    lib_ir,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_gnp_vintage,...
    lib_gnp_vintage,...
    lib_gnp,...
    lib_ea,...
    lib_cx,...
    lib_cx,...
    lib_fi,...
    lib_fi,...
    lib_fi,...
    lib_fi,...
    lib_fi,...
    lib_fi,...
    lib_fi,...
    lib_fi,...
    lib_ns_vintage,...
    lib_ns_vintage,...
    lib_ns_vintage,...
    lib_ns_vintage_gross,...
    lib_ns_vintage_gross,...
    lib_ns_vintage_gross,...
    lib_ns_vintage,...
    lib_ns_vintage,...
    lib_ns_vintage,...
    lib_ns_vintage_gross,...
    lib_ns_vintage_gross,...
    lib_ns_vintage_gross,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_tax,...
    lib_tax,...
    lib_tax,...
    lib_profits,...
    lib_profits,...
    lib_cc,...
    lib_cc,...
    lib_mip,...
    %'',...
    },...
    'src', {...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    %'',...
    },...
    'VariableNames', cell(1,51)...
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




if strcmp(freq, 'Q')
    %pass
elseif strcmp(freq, 'A')
    T = table2timetable(T);
    T = retime(T, 'yearly', to_Annual);
    T = timetable2table(T);
elseif strcmp(freq, 'M')
    T = table2timetable(T);
    T = retime(T, 'quarterly', 'mean');
    T = timetable2table(T);
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

%T = fillmissing(T, 'previous');


% Sample period
for period = sample_tab.period(end-2)'
%for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);
    
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

    T = outerjoin(T, EITB, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, UITB, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    
    EITB = T.EITB;
    UITB = T.UITB;
    T = removevars(T, {'EITB', 'UITB'});  

    if strcmp(freq, 'A')
        EITB = EITB*12;
        UITB = UITB*12;
    elseif strcmp(freq, 'Q')
        EITB = EITB*4;
        UITB = UITB*4;
    
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

    idx = ismember({data.series}, series_ns_fred);
    %NS = sum(T{:, [data(idx).VariableNames]}, 2); % current cost, net
    NS = T{:, [data(idx).VariableNames]}; % current cost, net; not summing, ipd...

    idx = ismember({data.series}, series_ns_fred_h);
    %NS_h = sum(T{:, [data(idx).VariableNames]}, 2); % real/historical cost, net
    NS_h = T{:, [data(idx).VariableNames]}; % real/historical, net; not summing, ipd...
        
    ipd_NS = NS./NS_h; % implicit price deflator
        
    idx = ismember({data.series}, series_inv);
    INV = sum(T{:, [data(idx).VariableNames]}, 2);
    
    
    % ***************************************
    % Gross capital stock
    
    idx = ismember({data.series}, strcat(series_ns_fred, '_gross'));
    %GS = sum(T{:, [data(idx).VariableNames]}, 2); % current cost, gross
    GS = sum(T{:, [data(idx).VariableNames]}, 2); % current cost, gross, net; not summing, ipd...

    idx = ismember({data.series}, strcat(series_ns_fred_h, '_gross'));
    %GS_h = sum(T{:, [data(idx).VariableNames]}, 2); % real/historical cost, gross
    GS_h = T{:, [data(idx).VariableNames]}; % real/historical cost, gross; not summing, ipd...
    
    
    % first difference in gross stock = gross investment flows
    dGS = [NaN(1, size(GS,2));...
        GS(2:end, :)-...
        GS(1:end-1, :)];

    dGS_h = [NaN(1, size(GS_h, 2));...
        GS_h(2:end, :)-...
        GS_h(1:end-1, :)];

    %[T(:, 'date'), array2table(dGS)]
    % TODO: procedure complete empty observations. This should not occur.
    % in the meantime this is just a check, so let it go.
    %CX = dGS; % need to use historical cost/real cost to get gross investment
    %CX = sum(dGS_h .* ipd_NS, 2); % summing
    
    % take 1.
    % result is not good. It was closer to Fama without the vintage.
    % I think gross investment needs some adjustment. not simply dGS
    % TODO: calc gross following the BEA document.
    % might need to get price deflator by asset type. :-(
    
    % take 2.
    % introduced implicit price deflator.
    % stdev is very close, average is off at levels
    % missing something in numerator or adding something not needed in
    % denominator. should be about twice...
    

    % ***************************************
    % Return on Capital, ROC
        
    idx = ismember({data.series}, series_pbt);
    pbt = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_mip);
    mip = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_tax);
    taxes = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_cc);
    CC = T{:, data(idx).VariableNames};    
    %CCAdj = pbt_tr_incl - pbt_tr - iva;
    
    % ***************************************
    % CRSP value weighted inedx

    idx = ismember({data.series}, series_crsp);
    CRSP = T{:, data(idx).VariableNames};
    
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
    
    
    ROC_num = pbt + mip + CC - taxes;
    ROC_denom = sum(NS, 2) + INV;
    
    
    ROC = ROC_num ./ ROC_denom;    
    CX_NS = CX ./ sum(NS, 2);
        
    dROC = [NaN;...
        ROC(2:end, :)-...
        ROC(1:end-1, :)];
    
    
    dCX_NS = [NaN;...
        CX_NS(2:end, :)-...
        CX_NS(1:end-1, :)];    

    dCRSP = [NaN;...
        CRSP(2:end, :)./...
        CRSP(1:end-1, :)];
    

   
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
    
    EITB = EITB/100*pcent_mult;
    UITB = UITB/100*pcent_mult;
    
    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;    


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

    EITB = EITB(idx_sample);
    UITB = UITB(idx_sample);
    
    clf
    subplot(1,2,1)
    hold on
    %plot(T.date(idx_sample), CX_NS(idx_sample))
    idx = ismember({data.series}, strcat(series_ns_fred, '_gross'));
    area(T.date(idx_sample), T{idx_sample, [data(idx).VariableNames]})
    plot(T.date(idx_sample), GS(idx_sample))
    T(idx_sample, [{'date'}, data(idx).VariableNames])
    yyaxis right
    plot(T.date(idx_sample), dGS(idx_sample))
    subplot(1,2,2)
    plot(T.date(idx_sample), CX_NS)
    
    sgtitle(sprintf('E[CX / NS]=%.5f; std[CX / NS]=%.5f',...
        mean(CX_NS, 'omitnan'), std(CX_NS, 'omitnan')))

    
    N = size(Pi, 1);
    %{
    [mean(CX_NS, 'omitnan'), std(CX_NS, 'omitnan');...
        mean(dCX_NS, 'omitnan'), std(dCX_NS, 'omitnan')]
    return
    %}
    tab4_ = [...
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
    tab4_ = array2table(tab4_, 'VariableNames', VariableNames);
    
    
    tab4_ = [cell2table(XNames, 'VariableNames', {'Variable (x)'}), tab4_];
    tab4_ = [array2table(repmat(max_date, size(tab4_, 1), 1), 'VariableNames', {'max_date'}), tab4_];
    tab4_ = [array2table(repmat(min_date, size(tab4_, 1), 1), 'VariableNames', {'min_date'}), tab4_];
    tab4 = [tab4; tab4_];
    

    
end


tab4
return
fpath_tab4 = fullfile(root, 'Fama1981', sprintf('tab4.%s.csv', freq));
writetable(tab4, fpath_tab4)

%{
UPDATE 2023_10_02:
Every single row in the table is different to Fama, 1981.

Fama uses the SCB 1977_01 and 1978_01 updates for GNP, or something
very close to that.

Fama uses log approximation.


I ran the RGNP with vintage series.
desc stats are close, not exact, but close.
Mean and standard deviation are almost exact.

TODO 2023_10_03: 
reform all data series to be vintage of roughly 1979_08
almost DONE! directories are marked vintage_yyyy_mm
within each section. Currently did net capital stock and gross
capital stock (the latter to calculate gross investment expenditure).

net and gross.
BOGZ1FL105013665A
BOGZ1FL105013265A
BOGZ1FL105012865A
This gives me CX_NS. I suggest checking the desc stats for that before
continuing on.

What's left?
0. vintage of real cost (constant dollar) capital stock, net & gross V
1. Inventories
2. pbt
3. irs figures
4. mip
5. cc allowance incl adj
6. tax

NS, net capital stock
series_ns_fred = {...
 'BOGZ1FL105013665A',... V
 'BOGZ1FL105013265A',... V
 'BOGZ1FL105012665A',... using BOGZ1FL105012865A = + BOGZ1FL105012665A + BOGZ1FL105012265A
 'BOGZ1FL105012265A'...  using BOGZ1FL105012865A = + BOGZ1FL105012665A + BOGZ1FL105012265A
};


CX, capital expenditures
uses ubpulished BEA data.
series_cx = {...
 'BOGZ1FA105013063A',... Nonresidential Structures; Annual transactions from BEA, Fixed Asset Table 4.7 Investment in Private Nonresidential Fixed Assets by Industry Group and Legal Form of Organization; line 39, Corporate nonfinancial, Structures.
 'BOGZ1FA105013023A',... Nonresidential Equipment; Annual transactions from BEA, Fixed Asset Table 4.7 Investment in Private Nonresidential Fixed Assets by Industry Group and Legal Form of Organization; line 38, Corporate nonfinancial, Equipment.
 'BOGZ1FA105012063A',... Residential Structures; Annual and quarterly transactions from unpublished BEA data
 'BOGZ1FA105012023A',... Residential Equipment; Annual and quarterly transactions from unpublished BEA data
};

PNFIA,  Private Nonresidential Fixed Investment
PRFI,   Private Residential Fixed Investment
No breakdown to legal form of incorporation (nonfinancial etc)

Can't use investment from NIPA tables. Must get changes in gross capital
stock in order to get gross investment; see perpetual inventory method:
``The perpetual inventory method cumulates past investment flows to
indirectly estimate the value of the stock.''

Depreciation and net capital stocks: Assets are carried in gross capital
stocks at their undepreciated value during the entire time they remain in
the stock. The value of these assets is depreciated to obtain net stocks,
which equal the difference between the cumulative value of gross investment
and cumulative depreciation.

TODO 2023_10_03:
calculate gross from nipa current series to have updated data.

UPDATE: 2023_10_04
so far I wan't able toreproduce the exact figures of CX/NS. My CX figure
are potentially missing a component accounting to roughly half of the
level. first differnces are closer with vintage data. I'm halting the
vintage project until I can safely compute CX. It seems from Fama's
footnote on the data that CX can be computedfrom tables 1-8 on Fixed
capital (first difference of gross times implicit price deflator etc), but
this has proved difficult. Waiting for advice.

TODO 2023_10_04:
1. get BLS vintage data following their email. FAIL, see email with BLS.
2. Summarize Fama results (print table 4 and 5).
3. Wait for advice on Israeli data to do the same for Israel.
4. Continue with other agenda: introduce the inflation factor returns in
regressions.



%}



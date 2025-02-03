%{

This script runs replicates 7 regressions in table 3, Fama (1981).

It does so with updated data, which is a bit different to the original.

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
lib_ea = fullfile(lib_data, 'economic_activity');  % economic activity
lib_cx = fullfile(lib_data, 'cx');  % capital expenditure
lib_tax = fullfile(lib_data, 'tax');  % capital expenditure
lib_ns = fullfile(lib_data, 'ns', 'level');  % capital expenditure
lib_cc = fullfile(lib_data, 'consumption_fixed_capital');  % capital consumption
lib_roc = fullfile(lib_data, 'roc');  % items for ROC_t
lib_mip = fullfile(lib_data, 'interest_paid');  % interest paid
lib_profits = fullfile(lib_data, 'profits');  % interest paid
lib_fi = fullfile(lib_data, 'fixed_investment');  % Fixed investment (gross)

lib_EITB = fullfile(root, 'EITB');
lib_UITB = fullfile(root, 'UITB');

lib_ir = fullfile(lib_data, 'ir');
lib_pop = fullfile(lib_data, 'population');
lib_consumption = fullfile(lib_data, 'consumption_usa');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
lag_inflation = 0;
%log_approx = true;
log_approx = false;
freq = 'Q';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
pcent_mult = 1; % set as 100 for percent

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
     datetime(1953,1,31), datetime(2007,7,31), datetime(1953,1,31),...
     datetime(1954,1,31), datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(2019,12,31), datetime(1977,12,31),...
    datetime(1976,12,31), datetime(1977,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});

%}

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

tab3_1 = [];
tab3_2 = [];
tab3_3 = [];
tab3_4 = [];
tab3_5 = [];
tab3_6 = [];
tab3_7 = [];


series_CPI = 'CPIAUCSL';
series_Rf = 'RF_FF';
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
    series_cx = {'BOGZ1FA105013063A', 'BOGZ1FA105013023A', 'BOGZ1FA105012063A', 'BOGZ1FA105012023A'};
    series_gnp = 'A791RX0Q048SBEA'; % Real gross national product per capita
    %series_gnp = 'GNPCA'; % Real Gross National Product

    % residential and nonresidential equipment and structures, capital stock
    series_ns_fred = {'BOGZ1FL105013665A', 'BOGZ1FL105013265A', 'BOGZ1FL105012665A', 'BOGZ1FL105012265A'}; % annual series
    series_roc_denom = [series_ns_fred, {'NCBICBA027N'}]; % series_ns_fred + inventories

elseif strcmp(freq, 'Q')
    series_pbt = 'A464RC1Q027SBEA'; % profit before tax
    series_tax = 'BOGZ1FL103178005Q'; % Nonfinancial Corporate Business; Total Taxes Payable; Liability, Level
    series_mip = 'BOGZ1FA106130001Q'; % Nonfinancial Corporate Business; Interest Paid, Transactions
    series_cc = 'BOGZ1FU106300015Q'; % Nonfinancial Corporate Business; Capital Consumption Allowance, Transactions
    series_inv = 'IABSNNCB'; % Nonfinancial Corporate Business; Inventories Excluding IVA, Current Cost Basis, Level
    series_cx_fred = 'BOGZ1FA105050005Q'; % Nonfinancial Corporate Business; Total Capital Expenditures, Transactions
    series_cx = {'BOGZ1FA105013063Q', 'BOGZ1FA105013023Q', 'BOGZ1FA105012063Q', 'BOGZ1FA105012023Q'};
    series_gnp = 'A791RX0Q048SBEA'; % Real gross national product per capita
    %series_gnp = 'GNPC96'; % Real Gross National Product
    
    % residential and nonresidential equipment and structures, capital stock
    series_ns_fred = {'RCSNNWMVBSNNCB', 'BOGZ1FL105013265Q', 'RCVSRNWMVBSNNCB', 'BOGZ1FL105012265Q'}; % quarterly series
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
    series_CPI,... Price index of all urban consumers
    series_Rf,...
    'BOGMBASE',...
    'CURRCIR',...
    'RESBALNS',...
    'GNPC96',... % real GNP
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
    'RCSNNWMVBSNNCB',... % Nonfinancial Corporate Business; Nonresidential Structures, Current Cost Basis, Level
    'BOGZ1FL105013265Q',... % Nonfinancial Corporate Business; Nonresidential Equipment, Current Cost Basis, Level
    'RCVSRNWMVBSNNCB',...  % Nonfinancial Corporate Business; Residential Structures, Current Cost Basis, Level
    'BOGZ1FL105012265Q',...  % Nonfinancial Corporate Business; Residential Equipment, Current Cost Basis, Level
    'BOGZ1FL105013665A',... % Nonfinancial Corporate Business; Nonresidential Structures, Current Cost Basis, Level
    'BOGZ1FL105013265A',... % Nonfinancial Corporate Business; Nonresidential Equipment, Current Cost Basis, Level
    'BOGZ1FL105012665A',...  % Nonfinancial Corporate Business; Residential Structures, Current Cost Basis, Level
    'BOGZ1FL105012265A',...  % Nonfinancial Corporate Business; Residential Equipment, Current Cost Basis, Level
    'IABSNNCB',... % Nonfinancial Corporate Business; Inventories Excluding IVA
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
    lib_ir,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_gnp,...
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
    %'',...
    },...
    'VariableNames', cell(1,35)...
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

T = fillmissing(T, 'previous');




% Sample period
for period = sample_tab.period(1)'

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
    NS = sum(T{:, [data(idx).VariableNames]}, 2);
    
    idx = ismember({data.series}, series_inv);
    INV = sum(T{:, [data(idx).VariableNames]}, 2);
    

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
    ROC_denom = NS + INV;
    
    
    ROC = ROC_num ./ ROC_denom;    
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
    
    EITB = EITB/100*pcent_mult;
    UITB = UITB/100*pcent_mult;

    % define leads and lags
    % T = sortrows(T, 'date', 'ascending');  % make sure T is ascending
    dPR_lag = [NaN; dPR(1:end-1)];
    dM_lag = [NaN; dM(1:end-1)];
    dCX_NS_lag = [NaN; dCX_NS(1:end-1)];
    dROC_lag = [NaN; dROC(1:end-1)];
    dRGNP_lag = [NaN; dRGNP(1:end-1)];
    Rf_lag = [NaN; Rf(1:end-1)];
    EITB_lag = [NaN; EITB(1:end-1)];
    UITB_lag = [NaN; UITB(1:end-1)];

    dPR_lead = [dPR(2:end); NaN];
    dM_lead = [dM(2:end); NaN];
    dCX_NS_lead = [dCX_NS(2:end); NaN];
    dROC_lead = [dROC(2:end); NaN];
    dRGNP_lead = [dRGNP(2:end); NaN];
    Rf_lead = [Rf(2:end); NaN];
    
    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

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

    EITB = EITB(idx_sample);
    EITB_lag = EITB_lag(idx_sample);
    UITB = UITB(idx_sample);
    UITB_lag = UITB_lag(idx_sample);
    
    N = size(Pi, 1);
    
    
    % **********************************
    % REGRESSION (1)
    Y = dCX_NS;
    X = [ones(N, 1), dROC, dROC_lag];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    tab3_1_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:5])];
    tab3_1 = [tab3_1; tab3_1_];

    % **********************************
    % REGRESSION (2)
    Y = dCX_NS;
    X = [ones(N, 1), dPR, dPR_lag];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    tab3_2_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:5])];
    tab3_2 = [tab3_2; tab3_2_];

    
    % **********************************
    % REGRESSION (3)
    Y = dROC;
    X = [ones(N, 1), dPR, dPR_lag];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    tab3_3_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:5])];
    tab3_3 = [tab3_3; tab3_3_];

    
    % **********************************
    % REGRESSION (4)
    Y = dCX_NS;
    X = [ones(N, 1), dROC, dROC_lag, dPR, dPR_lag];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    tab3_4_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:5])];
    tab3_4 = [tab3_4; tab3_4_];
    

    % **********************************
    % REGRESSION (5)
    Y = dCX_NS;
    X = [ones(N, 1), dROC, dPR_lag];
    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    tab3_5_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:5])];
    tab3_5 = [tab3_5; tab3_5_];
    

    % **********************************
    % REGRESSION (6)
    Y = dCX_NS;
    X = [ones(N, 1), dROC, dPR_lag, EITB, UITB, EITB_lag, UITB_lag];
    K = size(X, 2); % # of regressors
    
    idx_notna = all(~isnan([Y, X]), 2);
    Y = Y(idx_notna);
    X = X(idx_notna, :);
    N_ = numel(Y);
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N_-1)/(N_-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    tab3_6_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:5])];
    tab3_6 = [tab3_6; tab3_6_];

    
    % **********************************
    % REGRESSION (7)
    Y = dROC;
    X = [ones(N, 1), dPR, dPR_lag, EITB, UITB, EITB_lag, UITB_lag];
    %X = [ones(N, 1), EITB, UITB, EITB_lag, UITB_lag];
    K = size(X, 2); % # of regressors
    
    idx_notna = all(~isnan([Y, X]), 2);
    Y = Y(idx_notna);
    X = X(idx_notna, :);
    N_ = numel(Y);

    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    %[EstCoeffCov,se,coeff] = hac(X,Y, 'intercept', false);
    R_sq_adj = 1 - (1-STATS(1))*(N-1)/(N-K-1);
    
    
    sigma_hat = R'*R/(N-K);
    Sigma_hat = inv(D)*sigma_hat;
    t_stat = B./sqrt(diag(Sigma_hat));
    % not using hac
    tab3_7_ = [N, B', t_stat', STATS(1), R_sq_adj, sqrt(sigma_hat), arrayfun(@(lag) autocorr_(R, lag), [1:5])];
    tab3_7 = [tab3_7; tab3_7_];
    

    


    
end


tab3_1 = array2table(tab3_1, 'VariableNames', {...
    'nobs', 'constant', 'DROC_t', 'DROC_{t-1}',...
    'tstat_constant', 'tstat_DROC_t', 'tstat_DROC_{t-1}',...
    'R_sq', 'R_sq-adj', 's(epsilon)',...
    'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});

tab3_2 = array2table(tab3_2, 'VariableNames', {...
    'nobs', 'constant', 'DPR_t', 'DPR_{t-1}',...
    'tstat_constant', 'tstat_DPR_t', 'tstat_DPR_{t-1}',...
    'R_sq', 'R_sq-adj', 's(epsilon)',...
    'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});

tab3_3 = array2table(tab3_3, 'VariableNames', {...
    'nobs', 'constant', 'DPR_t', 'DPR_{t-1}',...
    'tstat_constant', 'tstat_DPR_t', 'tstat_DPR_{t-1}',...
    'R_sq', 'R_sq-adj', 's(epsilon)',...
    'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});

tab3_4 = array2table(tab3_4, 'VariableNames', {...
    'nobs', 'constant', 'DROC_t', 'DROC_{t-1}', 'DPR_t', 'DPR_{t-1}',...
    'tstat_constant', 'tstat_DROC_t', 'tstat_DROC_{t-1}', 'tstat_DPR_t', 'tstat_DPR_{t-1}',...
    'R_sq', 'R_sq-adj', 's(epsilon)',...
    'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});

tab3_5 = array2table(tab3_5, 'VariableNames', {...
    'nobs', 'constant', 'DROC_t', 'DPR_{t-1}',...
    'tstat_constant', 'tstat_DROC_t', 'tstat_DPR_{t-1}',...
    'R_sq', 'R_sq-adj', 's(epsilon)',...
    'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});

tab3_6 = array2table(tab3_6, 'VariableNames', {...
    'nobs', 'constant', 'DROC_t', 'DPR_{t-1}',...
    'EITB_{t-1}', 'UITB_t', 'EITB_{t-2}', 'UITB_{t-1}',...
    'tstat_constant', 'tstat_DROC_t', 'tstat_DPR_{t-1}',...
    'tstat_EITB_{t-1}', 'tstat_UITB_t', 'tstat_EITB_{t-2}', 'tstat_UITB_{t-1}',...
    'R_sq', 'R_sq-adj', 's(epsilon)',...
    'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});

tab3_7 = array2table(tab3_7, 'VariableNames', {...
    'nobs', 'constant', 'DPR_t', 'DPR_{t-1}',...
    'EITB_{t-1}', 'UITB_t', 'EITB_{t-2}', 'UITB_{t-1}',...
    'tstat_constant', 'tstat_DPR_t', 'tstat_DPR_{t-1}',...
    'tstat_EITB_{t-1}', 'tstat_UITB_t', 'tstat_EITB_{t-2}', 'tstat_UITB_{t-1}',...
    'R_sq', 'R_sq-adj', 's(epsilon)',...
    'rho_1', 'rho_2', 'rho_3', 'rho_4', 'rho_5'});


for tab_num = 1:7
    
fpath_tab3 = fullfile(root, 'Fama1981', sprintf('tab3_%d.%s.csv', tab_num, freq));
eval(sprintf('tab = tab3_%d;', tab_num))
tab = [sample_tab, tab];
writetable(tab, fpath_tab3)    
    
end





%%
clc


ROC(ROC==0) = NaN;
CX_NS(CX_NS==0) = NaN;
EITB

pos = [454   251   905   599];

size(dROC)
size(T)
clf
%idx = ~isnan(ROC);
%idx = ~isnan(CX_NS);
idx = ~isnan(EITB);

%plot(T.date(idx), ROC(idx), 'LineWidth', 1.5)
%plot(T.date(idx_sample), dROC)
%plot(T.date(idx), CX_NS(idx), 'LineWidth', 1.5)
%plot(T.date(idx_sample), dCX_NS)
plot(T.date(idx_sample), EITB)


yyaxis right
plot(T.date(T.date.Year>1951), PR(T.date.Year>1951), 'LineWidth', 1.5)
%plot(T.date(idx), PR(idx), 'LineWidth', 1.5)
%plot(T.date(idx_sample), dPR)

set(gcf, 'Color', 'w')
set(gca, 'FontSize', 12)
set(gcf, 'Position', pos)

fpath_fig = fullfile(root, 'gfx', 'EITB_usa.png');
%export_fig(fpath_fig);


mean(CX_NS(idx_sample), 'omitnan')






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
lib_gnp = fullfile(lib_data, 'gnp');  % gross national product
lib_ea = fullfile(lib_data, 'economic_activity');  % economic activity
lib_cx = fullfile(lib_data, 'cx');  % capital expenditure
lib_tax = fullfile(lib_data, 'tax');  % capital expenditure
lib_ns = fullfile(lib_data, 'ns', 'level');  % capital expenditure
lib_roc = fullfile(lib_data, 'roc');  % items for ROC_t
lib_ip = fullfile(lib_data, 'interest_paid');  % interest paid
lib_profits = fullfile(lib_data, 'profits');  % interest paid

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
freq = 'A';

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

tab2 = [];


series_CPI = 'CPIAUCSL';
series_Rf = 'RF_FF';
%series_ms = 'BOGMBASE';
series_ms = {'CURRCIR', 'RESBALNS'};
series_ea = 'INDPRO';

series_tax = 'BOGZ1FL143178005Q';
series_mip = 'BOGZ1FA106130001Q'; % monetary interest paid
series_pbt = 'A464RC1Q027SBEA'; % profit before tax
series_pbt_tr = 'BOGZ1FU106060005Q'; % profit before tax
series_pbt_tr_incl = 'BOGZ1FU106060035Q'; % profit before tax, incl iva and ccadj
series_iva = 'NCBIVDQ027S';

% note: there are discrepancies between NIPA and FRED. use NIPA
series_cx_fred = 'BOGZ1FA105050005Q';
series_cx_nipa = 'BOGZ1FA105050005A_NIPA';

%series_gnp = 'GNPC96';
series_gnp = 'A791RX0Q048SBEA';

% residential and nonresidential equipment and structures, capital stock
series_ns_fred = {'RCSNNWMVBSNNCB', 'BOGZ1FL105013265Q', 'RCVSRNWMVBSNNCB', 'BOGZ1FL105012265Q'};
series_roc_denom = [series_ns_fred, {'IABSNNCB'}]; % series_ns_fred + inventories
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
    'ESABSNNCB',... % Nonfinancial Corporate Business; Equipment
    'RCSNNWMVBSNNCB',... % Nonfinancial Corporate Business; Nonresidential Structures
    'RCVSRNWMVBSNNCB',...  % Nonfinancial Corporate Business; Residential Structures
    'BOGZ1FL105013265Q',...  % Nonfinancial Corporate Business; Residential Equipment
    'IABSNNCB',... % Nonfinancial Corporate Business; Inventories Excluding IVA
    'A464RC1Q027SBEA',... % Nonfinancial corporate business: Profits before tax (without IVA and CCAdj)
    'A466RC1A027NBEA',... % Nonfinancial corporate business: Profits after tax (without IVA and CCAdj)
    'BOGZ1FL143178005Q',... % Nonfinancial corporate business: Total tax liability
    'B470RC1Q027SBEA',... % Nonfinancial corporate business: Capital consumption adjustment
    'BOGZ1FA146060035Q',... % Nonfinancial Business; Corporate Profits Before Tax Including IVA and CCAdj, Transactions (SA)
    'BOGZ1FU106060005Q',... % Nonfinancial Corporate Business; Corporate Profits Before Tax *Excluding IVA and CCAdj,Transactions (NSA)
    'BOGZ1FU106060035Q',... % Nonfinancial Corporate Business; Corporate Profits Before Tax *Including IVA and CCAdj,Transactions (NSA)
    'NCBIVDQ027S',... % Nonfinancial Corporate Business; Inventory Valuation Adjustment (IVA),Transactions
    'BOGZ1FA106130001Q',... % Nonfinancial Corporate Business; Interest Paid, Transactions
    %'',...
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
    lib_profits,...
    lib_profits,...
    lib_profits,...
    lib_ip,...
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
    %'',...
    },...
    'VariableNames', cell(1,24)...
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

if strcmp(freq, 'M')
    %pass
elseif strcmp(freq, 'Q')
    T = table2timetable(T);
    T = retime(T, 'quarterly', 'lastvalue');
    T = timetable2table(T);
elseif strcmp(freq, 'A')
    T = table2timetable(T);
    T = retime(T, 'yearly', 'lastvalue');
    T = timetable2table(T);
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

T = fillmissing(T, 'previous');




% Sample period
for period = sample_tab.period'

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
    
    
    % ***************************************
    % Capital expenditures

    idx = ismember({data.series}, series_cx_nipa);
    CX = T{:, data(idx).VariableNames};

    % ***************************************
    % Net capital stock

    idx = ismember({data.series}, series_ns_fred);
    NS = sum(T{:, [data(idx).VariableNames]}, 2);

    % ***************************************
    % Return on Capital, ROC
        
    idx = ismember({data.series}, series_pbt);
    pbt = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_pbt_tr);
    pbt_tr = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_pbt_tr_incl);
    pbt_tr_incl = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_mip);
    mip = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_tax);
    taxes = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_iva);
    iva = T{:, data(idx).VariableNames};
    
    CCAdj = pbt_tr_incl - pbt_tr - iva;

    ROC_num = pbt + mip + CCAdj - taxes;
    ROC_denom = NS + iva;
    
    ROC = ROC_num ./ ROC_denom;
    
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

        Pi = log(Pi)*100;
        dPR = log(dPR)*100;
        dRGNP = log(dRGNP)*100;
        dM = log(dM)*100;
        dCX_NS = log(dCX_NS)*100;
        dROC = log(dROC)*100;

    else

        Pi = (Pi - 1)*100;
        dPR = (dPR - 1)*100;
        dRGNP = (dRGNP- 1)*100;
        dM = (dM - 1)*100;
        dCX_NS = (dCX_NS - 1)*100;
        dROC = (dROC- 1)*100;

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

    EITB = T.EITB(idx_sample);
    UITB = T.UITB(idx_sample);
    T = removevars(T, {'EITB', 'UITB'});  
    
    N = size(Pi, 1);
    tab2_a = ...
    [...
        corr(Pi, [dRGNP_lag, dRGNP, dRGNP_lead], 'rows', 'complete');...
        corr(Pi, [dPR_lag, dPR, dPR_lead], 'rows', 'complete');...
        corr(Pi, [dCX_NS_lag, dCX_NS, dCX_NS_lead], 'rows', 'complete');...
        corr(Pi, [dROC_lag, dROC, dROC_lead], 'rows', 'complete')...
    ];
    tab2_b = ...
    [...
        corr(EITB, [dRGNP_lag, dRGNP, dRGNP_lead], 'rows', 'complete');...
        corr(EITB, [dPR_lag, dPR, dPR_lead], 'rows', 'complete');...
        corr(EITB, [dCX_NS_lag, dCX_NS, dCX_NS_lead], 'rows', 'complete');...
        corr(EITB, [dROC_lag, dROC, dROC_lead], 'rows', 'complete')...
    ];
    tab2_c = ...
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
    tab2_a = [...
        table(min_dates, max_dates, 'VariableNames', {'min_date', 'max_date'}),...
        cell2table(Xnames, 'VariableNames', {'X_t'}),...
        cell2table(Ynames, 'VariableNames', {'Y_t+tau'}),...
        array2table(tab2_a, 'VariableNames', VariableNames)...
        ];
    
    Xnames = {'EITB', 'EITB', 'EITB', 'EITB'}';
    tab2_b = [...
        table(min_dates, max_dates, 'VariableNames', {'min_date', 'max_date'}),...
        cell2table(Xnames, 'VariableNames', {'X_t'}),...
        cell2table(Ynames, 'VariableNames', {'Y_t+tau'}),...
        array2table(tab2_b, 'VariableNames', VariableNames)...
        ];
    
    Xnames = {'BG', 'BG', 'BG', 'BG'}';
    tab2_c = [...
        table(min_dates, max_dates, 'VariableNames', {'min_date', 'max_date'}),...
        cell2table(Xnames, 'VariableNames', {'X_t'}),...
        cell2table(Ynames, 'VariableNames', {'Y_t+tau'}),...
        array2table(tab2_c, 'VariableNames', VariableNames)...
        ];
    

    tab2_ = [tab2_a; tab2_b; tab2_c];
    tab2 = [tab2; tab2_];
    
end



fpath_tab2 = fullfile(root, 'Fama1981', sprintf('tab2.%s.csv', freq));

if exist(fpath_tab2, 'file') == 2 && false
    tab2_ = readtable(fpath_tab2);
    tab2_.Properties.VariableNames = tab2.Properties.VariableNames;
    tab2 = [tab2; tab2_];
end
%}
%{
% save each iteration alone

fpath_tab2 = fullfile(root, 'Fama1981', 'tab2.csv');
counter = 0;
while exist(fpath_tab2, 'file') == 2
    fpath_tab1 = fullfile(root, 'Fama1981', sprintf('tab2_%d.csv', counter));
    counter = counter + 1;
end
        
%}
%tab1 = sortrows(tab1, {'min_date', 'max_date'});
writetable(tab2, fpath_tab2)




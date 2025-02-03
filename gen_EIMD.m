
%{

Table 1 in Fama (1981) specifies two types of regressions:
a. pp. 549
(6) Pi_t = -b_0 - b_1 dA_t - b_2 dRf_t + b_3 dM_t + eta_t (1953-77 period)
b. pp. 547
(4) Pi_t = alpha_{t-1} + TB_{t-1} + eta_t (1954-76 period)

In this script I ran the regression (6), where Pi_t is the inflation rate
for month, quarter or year t, and TB_{t-1} is the treasury bill rate
observed at the beginning of the month or quarter.

M_t and dA_t are annualized growth rates of the money base and 
real activity (e.g. industrial production, series_ea = 'INDPRO').

I get quite different R^2, so I'm guessing I'm not following Fama close
enough.



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
lib_ea = fullfile(lib_data, 'economic_activity');  % economic activity

lib_ir = fullfile(lib_data, 'ir');
lib_pop = fullfile(lib_data, 'population');
lib_consumption = fullfile(lib_data, 'consumption_usa');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


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
     datetime(1953,1,31), datetime(1953,1,31),...
     datetime(1954,1,31), datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(1977,12,31),...
    datetime(1976,12,31), datetime(1977,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

table1a = [];


series_CPI = 'CPIAUCSL';
if strcmp(freq, 'M')
    series_Rf = 'RF_FF'; % 1 month TBill, from Kenneth French website, Monthly, not annuazlied
elseif strcmp(freq, 'Q')
    series_Rf = 'TB3MS'; % 3-Month Treasury Bill Secondary Market Rate, Discount Basis, Percent, Monthly, Not Seasonally Adjusted
    lib_ir = fullfile(lib_ir, 'TBs');
elseif strcmp(freq, 'A')
    series_Rf = 'TB1YR'; % 1-Year Treasury Bill Secondary Market Rate, Discount Basis, Percent, Monthly, Not Seasonally Adjusted
    lib_ir = fullfile(lib_ir, 'TBs');
else
    error('freq must be M, Q or A')
end

%series_ms = 'BOGMBASE'; % first obs is Dec 1958...
series_ms = {'CURRCIR', 'RESBALNS'};
series_ea = 'INDPRO';
series_C = {'PCENDG', 'PCES'};
series_PD = {'PCEPINDG', 'PCEPIS'};


% src: BI, Bank of Israel old website; FAME
data = struct(...
    'series', {...
    series_CPI,... Price index of all urban consumers
    series_Rf,...
    'BOGMBASE',...
    'CURRCIR',...
    'RESBALNS',...
    series_ea,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_ir,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ea,...
    %'',...
    },...
    'src', {...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    %'',...
    },...
    'VariableNames', cell(1,6)...
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


    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames};

    idx = ismember({data.series}, series_ea);
    A = T{:, data(idx).VariableNames};
    dA = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    dA3 = [NaN(3, 1);...
        T{4:end,   data(idx).VariableNames}./...
        T{1:end-3, data(idx).VariableNames}];
    dA4 = [NaN(4, 1);...
        T{5:end,   data(idx).VariableNames}./...
        T{1:end-4, data(idx).VariableNames}];
    dA12 = [NaN(12,1);...
        T{13:end,   data(idx).VariableNames}./...
        T{1:end-12, data(idx).VariableNames}];

    dA_f = [dA(2:end); NaN(1, 1)];
    dA3_f = [dA3(4:end); NaN(3, 1)];
    dA4_f = [dA4(5:end); NaN(4, 1)];
    dA12_f = [dA12(13:end); NaN(12, 1)];
        

    
    idx = ismember({data.series}, series_ms);
    M = sum(T{:, [data(idx).VariableNames]}, 2);
    dM = [NaN;...
        M(2:end)./...
        M(1:end-1)];
    dM3 = [NaN(3, 1);...
        M(4:end)./...
        M(1:end-3)];
    dM4 = [NaN(4, 1);...
        M(5:end)./...
        M(1:end-4)];
    dM12 = [NaN(12,1);...
        M(13:end)./...
        M(1:end-12)];
    
    
    idx = ismember({data.series}, series_CPI);
    CPI = T{:, data(idx).VariableNames};
    Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    Pi3 = [NaN(3, 1);...
        T{4:end,   data(idx).VariableNames}./...
        T{1:end-3, data(idx).VariableNames}];
    Pi12 = [NaN(12, 1);...
        T{13:end,   data(idx).VariableNames}./...
        T{1:end-12, data(idx).VariableNames}];
    
    
    if log_approx

        Pi = log(Pi)*in_pcent;
        Pi3 = log(Pi3)*in_pcent;
        Pi12 = log(Pi12)*in_pcent;

        dA = log(dA)*in_pcent;
        dA3 = log(dA3)*in_pcent;
        dA12 = log(dA12)*in_pcent;
        dA_f = log(dA_f)*in_pcent;
        dA3_f = log(dA3_f)*in_pcent;
        dA12_f = log(dA12_f)*in_pcent;

        dM = log(dM)*in_pcent;
        dM3 = log(dM3)*in_pcent;
        dM12 = log(dM12)*in_pcent;

    else

        Pi = (Pi -1)*in_pcent;
        Pi3 = (Pi3 -1)*in_pcent;
        Pi12 = (Pi12 -1)*in_pcent;

        dA = (dA -1)*in_pcent;
        dA3 = (dA3 -1)*in_pcent;
        dA12 = (dA12 -1)*in_pcent;
        dA_f = (dA_f -1)*in_pcent;
        dA3_f = (dA3_f -1)*in_pcent;
        dA12_f = (dA12_f -1)*in_pcent;

        dM = (dM -1)*in_pcent;
        dM3 = (dM3 -1)*in_pcent;
        dM12 = (dM12 -1)*in_pcent;

    end



    CPI = CPI(idx_sample);
    Pi = Pi(idx_sample);
    Rf = Rf(idx_sample);
    dM = dM(idx_sample);
    dM4 = dM4(idx_sample);
    dM12 = dM12(idx_sample);
    dA = dA(idx_sample);
    dA4 = dA4(idx_sample);
    dA_f = dA_f(idx_sample);
    dA4_f = dA4_f(idx_sample);
    dA12 = dA12(idx_sample);
    dA12_f = dA12_f(idx_sample);
    
    
    N = size(Pi, 1);
        
    % estimate equation (6)
    Y = Pi;
    if strcmp(freq, 'M')
        X = [ones(N, 1), dM12, dA12, dA12_f];    
    elseif strcmp(freq, 'Q')
        X = [ones(N, 1), dM4, dA4, dA4_f];
    elseif strcmp(freq, 'A')
        X = [ones(N, 1), dM, dA, dA_f];
    end


    K = size(X, 2); % # of regressors
    
    D = X'*X; % inverse information matrix
    
    [B,BINT,R,~,STATS] = regress(Y, X);
    EIMD = X*B;
    UIMD = R;

    clf
    hold on
    plot(T.date(idx_sample), Pi)
    plot(T.date(idx_sample), EIMD)

    drawnow

    EIMD = [T(idx_sample, 'date'),...
        array2table(EIMD, 'VariableNames', {'EIMD'})];
    UIMD = [T(idx_sample, 'date'),...
        array2table(UIMD, 'VariableNames', {'UIMD'})];
    
    fname = sprintf('EIMD__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB = fullfile(root, 'EIMD', fname);
    writetable(EIMD, fpath_EITB)
    
    fname = sprintf('UIMD__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB = fullfile(root, 'UIMD', fname);
    writetable(UIMD, fpath_UITB)


    
end

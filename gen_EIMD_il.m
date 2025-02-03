
%{
Data for Israel.
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
lib_data = fullfile('/media', username, 'D', 'data', 'BOI');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')

lib_israel = fullfile(lib_data, 'Israel');

lib_cpi = fullfile(lib_israel);
lib_ms = fullfile(lib_israel, 'monagg');
lib_mbase = fullfile(lib_israel, 'money_base');
lib_reserves = fullfile(lib_israel, 'reserves');

lib_gnp = fullfile(lib_israel, 'GNP');  % gross national product

%lib_ea = fullfile(lib_israel, 'indprod');  % economic activity
lib_ea = fullfile(lib_israel);  % economic activity, monthly series

lib_ir = fullfile(lib_israel);
%lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M');
lib_pop = fullfile(lib_israel, 'POP');
lib_consumption = fullfile(lib_israel, 'C');
%lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
lag_inflation = 1;
in_pcent = 1;
log_approx = false;
pcent = 1; % set as 100 if want percent
freq = 'M'; % frequency to be used


series_CPI = 'CP_SA.M';
series_ms = 'ZA108.M'; % waiting for Elad from MOS
%series_ms = {'A346.M_E', 'ZA215.M'}; % ask Elad about RESNALNS, 
%series_ea = 'CLS11.TPR_C.Q_SA_CHAINED'; % excl diamonds
series_ea = 'mdd13.level.m'; % excl diamonds

% total investment, incl financial sector
series_cx = 'AM_INV_TOTALT_Q_N';

%series_gnp = 'GNP.Q_N';
series_gnp = 'GNP.Q_FP';

series_ns = {'AM_NSTOCK_TOTALT_Q_N'};
series_C = {'C_@DUR.Q_FP_SA'}; % excl durable, should contain services

% not needed
%series_PD = {'PCEPINDG', 'PCEPIS'}; % I have the real series, don't really
%need this. (implicit deflator)
%series_pop = 'DEM_POP_AVE_CHN.Q';



sample_tab = table(...
    [1:3]',...
    [datetime(1996,3,31), datetime(1996,3,31), datetime(2009,1,31),...
     ]',...
    [datetime(2019,12,31), datetime(2008,12,31), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});




%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});





if strcmp(freq, 'M')
    series_Rf = 'TSB_BAGR_MAKAM_01M.M'; % 1 month MAKAM, annualized
elseif strcmp(freq, 'Q')
    series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized
elseif strcmp(freq, 'A')
    series_Rf = 'TSB_BAGR_MAKAM_12M.M'; % 12 month MAKAM, annualized
else
    error('freq must be M, Q or A')
end



% src: BI, Bank of Israel old website; FAME
data = struct(...
    'series', {...
    series_CPI,... Price index of all urban consumers
    series_Rf,...
    'ZA108.M',...
    'A346.M_E',...
    'ZA215.M',...
    series_ea,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_ir,...
    lib_mbase,...
    lib_mbase,...
    lib_reserves,...
    lib_ea,...
    %'',...
    },...
    'src', {...
    'FAME',...
    'FAME',...
    'NONE',...
    'NONE',...
    'NONE',...
    'FAME',...
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
    fpath_EITB = fullfile(root, 'EIMD_il', fname);
    writetable(EIMD, fpath_EITB)
    
    fname = sprintf('UIMD__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB = fullfile(root, 'UIMD_il', fname);
    writetable(UIMD, fpath_UITB)



    
end

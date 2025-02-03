%{

In-sample procedure for TB model, FGLS estimation.
Data for Israel.

2024-06-11:
I've changed the estimation to a joint estimation of (theta, beta) via
Maximum likelihood (ML). For reference, I've kept computation of beta_FGLS
with the same method, but the ML is better.

The resulting beta is closer to -1 than FGLS. It also yields a low first
order autocorrelation for the unexpected inflation (UITB).

That it is still negative at the monthly frequency is beyond me.
This may mean that the changes in the short yields are driven by stronger
forces such as risk premia and/or liquidity factors.

it looks like all predictability is coming from change in the
expected real rate

2024-06-10:
I've added series_Rf, using the MAKAM yield curve. (MAKAM_yields)

Many times I get a newgative coefficient, which makes me think somehting is
wrong in this regression. theta does seem correct, and the resulting R_sq
etc shows this is something with the data itself.

Timing: I want to move to daily data, from the 16th to the 15th.
I need to make sure the timing i correct.

For example, MAKAM_yields15 is a monthly series, using the first business
day after the 15th of each month (say the 16th of Apr, 2024).
The inflation series dated 30th Apr regards Inflation between 1-30 Apr
and is announced on the 15th of May.

In order to match these I need to match these I need to push the MAKAM
yield to the end of April, which means not lagging inflation.

Lagging inflation adds one calendar month to inflation, so that Inflation
between 1-30 Apr, dated 30-Apr, will be dated 31-May.

Whatever I do, the prediction regression in Israel at monthly freq
produces a negative beta_FGLS. This has indeed has some prediction, but
the value of it is very low (I've made errors by including the t-residual
and effectively got exploded R_sq). The regression at annual freq
generates a coefficient of 1. The regression at quarterly freq
generates a coef of .25.

No matter if I use raw MAKAM, or the yield curve.

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
lib_israel = fullfile(lib_data, 'Israel');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')

lib_cpi = fullfile(lib_israel);
lib_ms = fullfile(lib_israel, 'monagg');
lib_mbase = fullfile(lib_israel, 'money_base');
lib_reserves = fullfile(lib_israel, 'reserves');

lib_gnp = fullfile(lib_israel, 'GNP');  % gross national product

lib_ea = fullfile(lib_israel, 'indprod');  % economic activity

lib_ir = fullfile(lib_israel);
%lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M');
lib_pop = fullfile(lib_israel, 'POP');
lib_consumption = fullfile(lib_israel, 'C');
%lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');

% flag to override period in downstream scripts: false by default
lag_inflation = 0;
log_approx = false;
pcent = 100; % set as 100 if want percent
% keep pcent as 100, otherwise computation limits on log-likelihood
freq = 'M'; % frequency to be used


series_CPI = 'CP_SA.M';
%series_ms = 'ZA108.M'; % waiting for Elad from MOS
series_ms = {'A346.M_E', 'ZA215.M'}; % ask Elad about RESNALNS, 
series_ea = 'CLS11.TPR_C_SA.Q_CHAINED'; % excl diamonds

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
    [1:6]',...
    [datetime(1996,3,31), datetime(1996,3,31), datetime(2009,1,31),...
    datetime(1996,3,31), datetime(2000,1,31), datetime(2000,1,31),...
     ]',...
    [datetime(2019,12,31), datetime(2008,12,31), datetime(2019,12,31),...
    datetime(2024,4,30), datetime(2024,4,30), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});



%series_Rf = 'MAKAM_yields';
series_Rfs = {'MAKAM_raw', 'MAKAM_raw_15', 'MAKAM_raw_15mean',...
        'MAKAM_yields', 'MAKAM_yields_15', 'MAKAM_yields_15mean',...
        'TELBOR', 'TELBOR_15', 'TELBOR_15mean',...
        };
for series_Rf = series_Rfs

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});
if strcmp(series_Rf, 'MAKAM_raw')
    % X-month MAKAM, annualized average of calendar month.
    % raw data gets the closest bond to X-months, which inserts noise, due
    % to liquidity as well as tick size. Not recommended.

    series_Rf = 'TSB_BAGR_MAKAM_01M.M'; % 1 month MAKAM, annualized
    %series_Rf = 'TSB_BAGR_MAKAM_02M.M'; % 2 month MAKAM, annualized
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized

    lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M', 'with_metadata');
    
elseif strcmp(series_Rf, 'MAKAM_raw_15')
    % X-month MAKAM, annualized from daily (not average) after 15 of the month
    % raw data gets the closest bond to X-months, which inserts noise, due
    % to liquidity as well as tick size. Not recommended.

    series_Rf = 'TSB_BAGR_MAKAM_01M.M_15'; % 1 month
    %series_Rf = 'TSB_BAGR_MAKAM_02M.M_15'; % 2 month
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M_15'; % 3 month

    lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM_15', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_raw_15mean')
    % X-month MAKAM, annualized from daily (average) after 15 of the month
    % to the 15th of the next month.
    % raw data gets the closest bond to X-months, which inserts noise, due
    % to liquidity as well as tick size. Not recommended.

    series_Rf = 'TSB_BAGR_MAKAM_01M.M_15mean'; % 1 month
    %series_Rf = 'TSB_BAGR_MAKAM_02M.M_15mean'; % 2 month
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M_15mean'; % 3 month

    lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM_15mean', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_yields')
    % X-month MAKAM, annualized average of calendar month --
    % yield curve estimation using NS model.

    series_Rf = 'MAKAM_yields_M01.M'; % 1 month
    lib_ir = fullfile(lib_israel, 'MAKAM_yields', 'M', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_yields_15')
    % X-month MAKAM, annualized from daily (not average) after 15 of the
    % month --
    % yield curve estimation using NS model.

    series_Rf = 'MAKAM_yields_M01.M_15'; % 1 month
    lib_ir = fullfile(lib_israel, 'MAKAM_yields_15', 'with_metadata');

elseif strcmp(series_Rf, 'MAKAM_yields_15mean')
    % X-month MAKAM, annualized from daily (average) after 15 of the month
    % to the 15th of the next month --
    % yield curve estimation using NS model.

    series_Rf = 'MAKAM_yields_M01.M_15mean'; % 1 month
    lib_ir = fullfile(lib_israel, 'MAKAM_yields_15mean', 'with_metadata');

elseif strcmp(series_Rf, 'TELBOR')
    % X-month MAKAM, annualized average of calendar month --
    % TEBOR rates

    series_Rf = 'BL.TELBOR_01M.M'; % 1 month
    lib_ir = fullfile(lib_israel, 'TELBOR', 'M', 'with_metadata');

elseif strcmp(series_Rf, 'TELBOR_15')
    % X-month MAKAM, annualized from daily (not average) after 15 of the
    % month --
    % yield curve estimation using NS model.

    series_Rf = 'BL.TELBOR_01M.M_15'; % 1 month
    lib_ir = fullfile(lib_israel, 'TELBOR_15', 'with_metadata');
    
elseif strcmp(series_Rf, 'TELBOR_15mean')
    % X-month MAKAM, annualized from daily (average) after 15 of the month
    % to the 15th of the next month --
    % yield curve estimation using NS model.

    series_Rf = 'BL.TELBOR_01M.M_15mean'; % 1 month
    lib_ir = fullfile(lib_israel, 'TELBOR_15mean', 'with_metadata');

else
    error('Unknown %s', series_Rf)
end

if strcmp(freq, 'M')
    % pass
elseif strcmp(freq, 'Q')
    series_Rf = strrep(series_Rf, '01', '03');
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M_15';
    %series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized
elseif strcmp(freq, 'A')
    series_Rf = strrep(series_Rf, '01', '12');
    %series_Rf = 'TSB_BAGR_MAKAM_12M.M'; % 12 month MAKAM, annualized
else
    error('freq must be M, Q or A')
end


% src: BI, Bank of Israel old website; FAME
data = struct(...
    'series', {...
    'CP_NSA.M',... CPI, not seasonally adjusted, BLS
    'CP_SA.M',... Price index of all urban consumers
    series_Rf,...
    'ZA108.M',...
    'A346.M_E',...
    'ZA215.M',...
    series_ea,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
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
    'FAME',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    %'',...
    },...
    'VariableNames', cell(1,7)...
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
    [T_, metadata_] = load_data_files(data(i).series,...
        lib_data_series, data(i).src);
    data(i).VariableNames = T_.Properties.VariableNames(2);

    if strcmp(data(i).series, series_Rf)
        T_.date.Day = 1;
        T_.date = T_.date+calmonths(1);
    end

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


%{
% error in BOI CP_SA series... in Jan-2000.
% seasonally adjust CP_NSA for period not available from CBS, adding 12
% months for sa-estimation
idx_CP_SA = 1:find(~isnan(T.CP_SA0x2EM), 1, 'first')+12;
CP_SA = sa_adj(T.CP_NSA0x2EM(idx_CP_SA), 12);
% set CP_SA as the official CBS series of the CPI
%T.CP_SA = T.CP_SA0x2EM; % error in the official BOI series
%}

idx_CP_SA = T.date.Year > 1990;
CP_SA = sa_adj(T.CP_NSA0x2EM(idx_CP_SA), 12);


T.CP_SA(idx_CP_SA) = CP_SA;

%{
% complete this series with the above CP_SA
T.CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1) = ...
    CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1);
%}



grid_h = 1:12;
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
    grid_h = 1:4;
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

%T = fillmissing(T, 'previous');



params = containers.Map();
params('constant') = true;
params('weight') = false;
params('freq') = freq;


% Sample period
for period = sample_tab.period(1:end)'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);

    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;

    idx = ismember({data.series}, series_Rf);
    T.Rf = T{:, data(idx).VariableNames}/100;
    if strcmp(freq, 'M')
        T.Rf = T.Rf/12; % data in annual terms -- to monthly
    elseif strcmp(freq, 'Q')
        T.Rf = T.Rf/4;
    elseif strcmp(freq, 'A')
        % pass
    else
        error('freq must be M, Q or A')
    end
    
    T.Rf_lag = lag(T.Rf, 1);
    T.Rf_lead = lead(T.Rf, 1);
    
    T.Pi = [NaN;...
        T.CP_SA(2:end, :)./...
        T.CP_SA(1:end-1, :)];

    %{
    idx = ismember({data.series}, series_CPI);
    T.Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    %}
    
    if log_approx
        T.Pi = log(T.Pi)*pcent;
    else
        T.Pi = (T.Pi -1)*pcent;
    end
    T.Rf = T.Rf*pcent;

    T_ = T(idx_sample, :);

    summary_stats = [...
        [arrayfun(@(h) autocorr_(T_.Pi, h), grid_h), mean(T_.Pi), std(T_.Pi)];...
        [arrayfun(@(h) autocorr_(T_.Rf, h), grid_h), mean(T_.Rf), std(T_.Rf)]...
        ];
    array2table(summary_stats, 'VariableNames',...
        [strcat('rho_', arrayfun(@num2str, grid_h, 'UniformOutput', false)),...
        {'mean', 'std'}])

    lib_load_EITB = fullfile(root, 'EITB_il', series_Rf);
    lib_load_UITB = fullfile(root, 'UITB_il', series_Rf);


    fname = sprintf('EITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_EITB = fullfile(lib_load_EITB, fname);
    if exist(fpath_EITB, 'file') ~= 2
        continue
    end

    EITB = readtable(fpath_EITB);
    
    fname = sprintf('UITB__%s-%s.%s.csv', min_date, max_date, freq);
    fpath_UITB = fullfile(lib_load_UITB, fname);
    UITB = readtable(fpath_UITB);


    clf
    hold on
    plot(UITB.date, UITB.UITB)
    
    N = size(UITB, 1);
    h=12;
    ma12 = arrayfun(@(t) mean(UITB.UITB(t:t+h)), 1:N-h)';
    plot(UITB.date(1+h:end), ma12)
    title(sprintf('E(UITB)=%.3f', mean(UITB.UITB)))
    ylabel('Percent')
    tit = sprintf('%s: %s -- %s', series_Rf, min_date, max_date);
    tit = str4fig(tit);
    title([tit,...
        'Maximum likelihood',...
        sprintf('Mean = %.3f; Stdev = %.3f; ACF(1)=%.3f',...
        mean(UITB.UITB, 'omitnan'),...
        std(UITB.UITB, 'omitnan'),...
        autocorr_(UITB.UITB, 1))])

    legend({'Unexpected inflation (UITB)', '12 months moving average (UITB)'})
    
    drawnow
    fname_fig = sprintf('%s -- %s.jpg', min_date, max_date);
    lib_fig = fullfile(root, 'gfx', 'IL', 'UITB', series_Rf);
    if exist(lib_fig, 'dir') ~= 2
        mkdir(lib_fig)
    end

    fpath_fig = fullfile(lib_fig, fname_fig);
    export_fig(fpath_fig)


    
end
clear metadata T_ T
end


return

%%

period = 1;
min_date = sample_tab.min_date(period);
max_date = sample_tab.max_date(period);


fpath = fullfile(matlab_dir, 'NRC', 'arima_101_il.csv')
arima_101 = readtable(fpath);
arima_101 = renamevars(arima_101, 'Pi_hat', 'Pi_hat_101');
arima_101.date = arima_101.date+1;
arima_101.date = arima_101.date - calmonths(1);
arima_101.date = arima_101.date-1;

EITB = readtable(fpath_EITB);
UITB = readtable(fpath_UITB);

TB_model = join(EITB, UITB);

common_dates = intersect(arima_101.date, TB_model.date);
idx_arima = ismember(arima_101.date, common_dates);
idx_TB = ismember(TB_model.date, common_dates);
idx_T = ismember(T.date, common_dates);

clf
hold on
%plot(T.date(idx_T), T.Pi(idx_T), 'r-')
plot(arima_101.date(idx_arima), arima_101.Pi_hat_101(idx_arima), 'k-')
plot(TB_model.date(idx_TB), TB_model.EITB(idx_TB), 'r--')

sqrt(mean((arima_101.epsilon_101(idx_arima)).^2))
sqrt(mean((TB_model.UITB(idx_TB)).^2))

%legend({'\pi', 'EIAR', 'EITB'})
legend({'EIAR', 'EITB'})
ylabel('Percent')


corr(TB_model.EITB(idx_TB), arima_101.Pi_hat_101(idx_arima))
corr(TB_model.UITB(idx_TB), arima_101.epsilon_101(idx_arima))





%{
Professional forecasters of Israeli inflation, non-seasonally adjusted.
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

li = 2; % forecast horizon in months


fpath = fullfile(lib_israel, 'UITB.xlsx');
UITB = readtable(fpath);

fpath = fullfile(matlab_dir, 'NRC', 'arima_101_il.csv')
arima_101 = readtable(fpath);
arima_101 = renamevars(arima_101, 'Pi_hat', 'Pi_hat_101');

fpath = fullfile(matlab_dir, 'NRC', 'arima_011_il.csv')
arima_011 = readtable(fpath);
arima_011 = renamevars(arima_011, 'Pi_hat', 'Pi_hat_011');


fname = 'CP_NSA.M.xlsx';
fpath = fullfile(lib_israel, fname);
CP_NSA = load_FAME(fpath);
CP_NSA.Pi_NSA = [NaN; CP_NSA{2:end, end}./CP_NSA{1:end-1, end}];
CP_NSA.Pi_NSA = (CP_NSA.Pi_NSA-1)*100;
%CP_NSA.Pi_NSA = log(CP_NSA.Pi_NSA);

fname = 'CP_SA.M.xlsx';
fpath = fullfile(lib_israel, fname);
CP_SA = load_FAME(fpath);
CP_SA.Pi_SA = [NaN; CP_SA{2:end, end}./CP_SA{1:end-1, end}];
CP_SA.Pi_SA = (CP_SA.Pi_SA-1)*100;


CP_SA = innerjoin(CP_SA, arima_101(:, {'date', 'epsilon_101'}));
CP_SA = innerjoin(CP_SA, arima_011(:, {'date', 'epsilon_011'}));

%{
% check that Pi is the same
%CP_SA = innerjoin(CP_SA, arima_101);

CP_NSA.date = CP_NSA.date-calmonths(1);
CP_SA.date = CP_SA.date-calmonths(1);
%}

CP_NSA = table2timetable(CP_NSA);
CP_SA = table2timetable(CP_SA);

CP_NSA = retime(CP_NSA, 'daily', 'previous');
CP_SA = retime(CP_SA, 'daily', 'previous');

CP_NSA = timetable2table(CP_NSA);
CP_SA = timetable2table(CP_SA);





haz_type = 'AVG';
if li == 2.5
    freq = 'D';
    series_type = '';
    haz_horizon = '02M';

    lib_HAZ_PI = fullfile(lib_israel, ['HAZ_PI', series_type], freq);
    fname_HAZ = ['HAZ_MAD_', haz_type, '_', haz_horizon, sprintf('.%s', freq), series_type, '.csv'];
elseif li == 2.0
    freq = 'M';
    haz_horizon = '01M';
    series_type = '_15';
    
    lib_HAZ_PI = fullfile(lib_israel, ['HAZ_PI', series_type]);
    fname_HAZ = ['HAZ_MAD_', haz_type, '_', haz_horizon, sprintf('.%s', freq), series_type, '.csv'];
elseif li == 1.5
    freq = 'D';
    series_type = '';
    haz_horizon = '01M';

    lib_HAZ_PI = fullfile(lib_israel, ['HAZ_PI', series_type], freq);
    fname_HAZ = ['HAZ_MAD_', haz_type, '_', haz_horizon, sprintf('.%s', freq), series_type, '.csv'];
elseif li == 1.0
    freq = 'M';
    haz_horizon = '00M';
    series_type = '_15';
    
    lib_HAZ_PI = fullfile(lib_israel, ['HAZ_PI', series_type]);
    fname_HAZ = ['HAZ_MAD_', haz_type, '_', haz_horizon, sprintf('.%s', freq), series_type, '.csv'];
elseif li == 0.5
    freq = 'D';
    series_type = '';
    haz_horizon = '00M';

    lib_HAZ_PI = fullfile(lib_israel, ['HAZ_PI', series_type], freq);
    fname_HAZ = ['HAZ_MAD_', haz_type, '_', haz_horizon, sprintf('.%s', freq), series_type, '.csv'];
end

horizon_months = str2double(haz_horizon(1:end-1));

fpath_HAZ = fullfile(lib_HAZ_PI, fname_HAZ);
T = readtable(fpath_HAZ);
haz_series = T.Properties.VariableNames{end};

if li == 2.5
    % compared to eom UITB
    T.date = T.date + 1 + calmonths(horizon_months);
    T = T(T.date.Day == 1, :);
    T.date = T.date - 1;
    series_AR_oos = 'UIAR3';
elseif li == 2.0
    T = T(2:end, :);
    T.date.Day = 1;
    T.date = T.date + calmonths(2);
    T.date = T.date - 1;
    series_AR_oos = 'UIAR2';
elseif li == 1.5
    % compared to eom UITB
    T.date = T.date + 1 + calmonths(horizon_months);
    T = T(T.date.Day == 1, :);
    T.date = T.date - 1;
    series_AR_oos = 'UIAR2';
elseif li == 1.0
    T = T(2:end, :);
    T.date.Day = 1;
    T.date = T.date + calmonths(1);
    T.date = T.date - 1;
    series_AR_oos = 'UIAR1';
elseif li == 0.5
    % compared to eom UITB
    T.date = T.date + 1 + calmonths(horizon_months);
    T = T(T.date.Day == 1, :);
    T.date = T.date - 1;
    series_AR_oos = 'UIAR1';
end




%{
% to end of month
T.date = T.date+calmonths(horizon_months);
T.date = eom(T.date);
%}


T = innerjoin(T, CP_NSA);
T = innerjoin(T, CP_SA);




err_NSA = T.Pi_NSA-T.(haz_series);
err_SA = T.Pi_SA-T.(haz_series);

T.UIPF = err_NSA;
T.EIPF = T.Pi_SA - err_NSA;


% consider the timing of UIPF, 1.0 month, or 1.5 months
% 1.0 month --> use series_type = '_15' & haz_horizon = '00M', freq = 'M'
% 1.5 months --> use series_type = '' & haz_horizon = '01M', freq = 'D'
EIPF = T(:, {'date', 'Pi_SA', 'EIPF'});
UIPF = T(:, {'date', 'Pi_SA', 'UIPF'});

EIPF.date = datestr(EIPF.date, 'yyyy-mm-dd');
UIPF.date = datestr(UIPF.date, 'yyyy-mm-dd');

fpath_EIPF = fullfile(lib_israel, sprintf('EIPF_%.1f.xlsx', li));
fpath_UIPF = fullfile(lib_israel, sprintf('UIPF_%.1f.xlsx', li));

writetable(EIPF, fpath_EIPF)
writetable(UIPF, fpath_UIPF)

if li == 1
    fpath_EIPF = fullfile(lib_israel, 'EIPF.xlsx');
    fpath_UIPF = fullfile(lib_israel, 'UIPF.xlsx');
    
    writetable(EIPF, fpath_EIPF)
    writetable(UIPF, fpath_UIPF)
    
end


return
% compare with UITB
T = innerjoin(T, UITB(:, {'date', 'UITB'}));


err_NSA = T.Pi_NSA-T.(haz_series);
err_SA = T.Pi_SA-T.(haz_series);

N = size(T,1);
clf
hold on
h=0;
plot(T.date(1+h:end), arrayfun(@(t) mean(err_NSA(t:t+h)), 1:N-h)')
h=24;
plot(T.date(1+h:end), arrayfun(@(t) mean(err_NSA(t:t+h)), 1:N-h)')

title("Professional Forecasters' Error")

clf
grid_h = 1:12;
acf_HAZ = arrayfun(@(h) autocorr_(T.(haz_series), h), grid_h);
acf_Pi = arrayfun(@(h) autocorr_(T.Pi_NSA, h), grid_h);
bar(grid_h, [acf_Pi; acf_HAZ])

legend({'\pi', 'PF'}, 'Location','north')
title('Correlogram')

clf
hold on
grid_h = 1:12;
acf_err_NSA = arrayfun(@(h) autocorr_(err_NSA, h), grid_h);
acf_err_SA = arrayfun(@(h) autocorr_(err_SA, h), grid_h);
bar(grid_h, [acf_err_NSA; acf_err_SA])

legend({'err-NSA', 'err-SA'}, 'Location','north')
title('Correlogram')




%plot(T.date, err_SA)

[mean(err_NSA), mean(T.epsilon_101), mean(T.UITB)]
[autocorr_(err_NSA, 1), autocorr_(T.epsilon_101, 1), autocorr_(T.UITB, 1)]


corrcoef([err_NSA, T.epsilon_101, T.UITB])
sqrt(mean([err_NSA, T.epsilon_101, T.UITB].^2)) % RMSE
%corrcoef([err_NSA, T.epsilon_101, T.epsilon_011, T.UITB])
%sqrt(mean([err_NSA, T.epsilon_101, T.epsilon_011, T.UITB].^2)) % RMSE

%}





%%











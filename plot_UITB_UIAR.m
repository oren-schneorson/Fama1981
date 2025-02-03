% plot the ma and acf of UITB and UIAR

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

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')

% ARMA params
p = 1;
d = 0;
q = 1;

forecast_lag = 0; % in-sample
%forecast_lag = 1; % out-of-sample, 1-step ahead
%forecast_lag = 2; % out-of-sample, 2-step ahead

%series_Rf = 'MAKAM_yields_M01.M_15';
series_Rf = 'MAKAM_yields_M01.M';
%series_Rf = 'BL.TELBOR_01M.M';
if strcmp(series_Rf, 'BL.TELBOR_01M.M') && forecast_lag > 0
    min_date = datetime(2001, 10, 31);
    max_date = datetime(2019, 12, 31);
elseif strcmp(series_Rf, 'MAKAM_yields_M01.M') && forecast_lag > 0
    min_date = datetime(2000, 1, 31);
    max_date = datetime(2019, 12, 31);
else
    min_date = datetime(2000, 1, 31);
    max_date = datetime(2019, 12, 31);
end


fname_UIAR = sprintf('arima_%d%d%d_il.csv',p,d, q);
fpath_UIAR = fullfile(root, '..', 'NRC', fname_UIAR);

T = readtable(fpath_UIAR);
idx = T.date >= min_date & T.date <= max_date;
T = T(idx, :);

if forecast_lag > 0
    UIAR = T{:, sprintf('UIAR%d_%d%d%d', forecast_lag, p, d, q)};
    EIAR = T{:, sprintf('EIAR%d_%d%d%d', forecast_lag, p, d, q)};
    Pi_AR = T.Pi;
else
    idx = T.date >= min_date & T.date <= max_date;
    UIAR = T{idx, sprintf('epsilon_%d%d%d', p, d, q)};
    EIAR = T{idx, 'Pi_hat'};
    Pi_AR = T.Pi(idx);
    dates_AR = T.date(idx);

end




if forecast_lag > 0
    fname_UITB = [series_Rf, '.xlsx'];
    %fname_UITB = 'MAKAM_yields_M01.M.xlsx';
    fpath_UITB = fullfile(root, 'ML_FGLS_il_ZRD_oos', 'Rf', fname_UITB);
    
    T = readtable(fpath_UITB);
    T.date = T.date + 1;
    T.date = T.date - calmonths(1);
    T.date = T.date - 1;
    
    idx = T.date >= min_date & T.date <= max_date;
    T = T(idx, :);
    
    UITB = T.err_t_1__oos;
    UITB1 = T.err_1_t_1__oos;
    
    EITB = T.EITB_t_oos;
    EITB1 = T.EITB_1_t_oos;
    Pi_TB = T.Pi_t_1_;

else
    if strcmp(series_Rf, 'MAKAM_yields_M01.M')
        series_min_date = datetime(1996, 3, 31);
        series_max_date = datetime(2024, 3, 31);
    elseif strcmp(series_Rf, 'BL.TELBOR_01M.M')
        series_min_date = datetime(1999, 11, 30);
        series_max_date = datetime(2024, 3, 31);
    end

    fname_EITB = sprintf('EITB__%s-%s.M.xlsx', series_min_date, series_max_date);
    fname_EITB1 = sprintf('EITB1__%s-%s.M.xlsx', series_min_date, series_max_date);
    fname_UITB = sprintf('UITB__%s-%s.M.xlsx', series_min_date, series_max_date);
    fname_UITB1 = sprintf('UITB1__%s-%s.M.xlsx', series_min_date, series_max_date);


    fpath_EITB = fullfile(root, 'EITB_il', series_Rf, fname_EITB);
    fpath_EITB1 = fullfile(root, 'EITB_il', series_Rf, fname_EITB1);
    fpath_UITB = fullfile(root, 'UITB_il', series_Rf, fname_UITB);
    fpath_UITB1 = fullfile(root, 'UITB_il', series_Rf, fname_UITB1);
    
    EITB = readtable(fpath_EITB);
    EITB.date = datetime(EITB.date);
    EITB1 = readtable(fpath_EITB1);
    EITB1.date = datetime(EITB1.date);
    UITB = readtable(fpath_UITB);
    UITB.date = datetime(UITB.date);
    UITB1 = readtable(fpath_UITB1);
    UITB1.date = datetime(UITB1.date);

    
    EITB.date = EITB.date + 1;
    EITB.date = EITB.date - calmonths(1);
    EITB.date = EITB.date - 1;

    EITB1.date = EITB1.date + 1;
    EITB1.date = EITB1.date - calmonths(1);
    EITB1.date = EITB1.date - 1;

    UITB.date = UITB.date + 1;
    UITB.date = UITB.date - calmonths(1);
    UITB.date = UITB.date - 1;

    UITB1.date = UITB1.date + 1;
    UITB1.date = UITB1.date - calmonths(1);
    UITB1.date = UITB1.date - 1;

    idx = EITB.date >= min_date & EITB.date <= max_date;
    EITB = EITB(idx, :);
    idx = EITB1.date >= min_date & EITB1.date <= max_date;
    EITB1 = EITB1(idx, :);

    idx = UITB.date >= min_date & UITB.date <= max_date;
    UITB = UITB(idx, :);
    idx = UITB1.date >= min_date & UITB1.date <= max_date;
    UITB1 = UITB1(idx, :);
    
    %[EITB, UITB(:, 'UITB'), array2table(EITB.EITB+UITB.UITB)]
    %return

    EITB = EITB.EITB;
    EITB1 = EITB1.EITB1;
    UITB = UITB.UITB;
    UITB1 = UITB1.UITB1;
    
    Pi_TB = UITB+EITB;
    if abs(sum((UITB+EITB)-(EITB1+UITB1))) > 1e-10
        error('TB and TB1 series not in sync.')
    end

end

if abs(sum(Pi_AR-Pi_TB)) > 1e-10
    [array2table(dates_AR), array2table([Pi_AR,Pi_TB])]
    error('Inflation series not in sync.')
end

rmse(Pi_AR, EIAR)
rmse(Pi_TB, EITB)
rmse(Pi_TB, EITB1)

corr([UITB1, UIAR])
return


N = size(T, 1);

clf
hold on
hs = [36];
Colors = {'b', 'r'};
LineStyles = {'-', '--', ':'};
LineWidths = [1, 1, 2];
for h = hs
    ma_UIAR = arrayfun(@(t) mean(UIAR(t:t+h)), 1:N-h)';
    ma_UITB = arrayfun(@(t) mean(UITB1(t:t+h)), 1:N-h)';
    Color = Colors{hs==h};
    LineStyle = LineStyles{hs==h};
    LineWidth = LineWidths(hs==h);
    plot(T.date(1+h:end), ma_UIAR,...
        'Color', Colors{1}, 'LineStyle', LineStyle, 'LineWidth', LineWidth)
    plot(T.date(1+h:end), ma_UITB,...
        'Color', Colors{2}, 'LineStyle', LineStyle, 'LineWidth', LineWidth)
    %autocorr_(ma_UIAR, 1)
    %autocorr_(ma_UITB, 1)
    std(ma_UIAR)
    std(ma_UIAR)/std(UIAR)

    std(ma_UITB)
    std(ma_UITB)/std(UITB1)
end

leg = {'UIAR', 'UITB'};
legend(leg, 'Location', 'north')

%title(sprintf('E(UITB)=%.3f', mean(UITB)))
%}

set(gcf, 'Color', 'w')
fpath_fig = fullfile(root, 'gfx', 'UIAR_UITB_il_ma.jpg');
export_fig(fpath_fig)

return
%%
clc; clf
hs = 1:12;
acf_UIAR = arrayfun(@(h) autocorr_(UIAR, h), hs);
acf_UITB = arrayfun(@(h) autocorr_(UITB1, h), hs);

sum(acf_UIAR)
sum(acf_UITB)

arrayfun(@(h) autocorr_(UIAR, h), hs)
arrayfun(@(h) autocorr_(UITB1, h), hs)
bar(hs, [acf_UIAR; acf_UITB])

set(gcf, 'Color', 'w')
legend({'UIAR', 'UITB'})

fpath_fig = fullfile(root, 'gfx', 'UITB_UIAR_il_acf.jpg');
export_fig(fpath_fig)


%%
clf
hold on
plot(T.date, Pi_TB, 'k--', 'LineWidth',1.5)
plot(T.date, EIAR, 'b')
plot(T.date, EITB1, 'r')

legend({'\pi', 'EIAR', 'EITB1'})
set(gcf, 'Color', 'w')


fpath_fig = fullfile(root, 'gfx', 'pi_EIAR_EITB.jpg');
export_fig(fpath_fig)

%%
clf
hold on
plot(T.date, Pi_TB, 'k--', 'LineWidth',1.5)
plot(T.date, EIAR, 'b')

legend({'\pi', 'EIAR', })
set(gcf, 'Color', 'w')


fpath_fig = fullfile(root, 'gfx', 'pi_EIAR.jpg');
export_fig(fpath_fig)


%%
clf
hold on
plot(T.date, Pi_TB, 'k--', 'LineWidth',1.5)
plot(T.date, EITB1, 'r')

legend({'\pi', 'EITB', })
set(gcf, 'Color', 'w')


fpath_fig = fullfile(root, 'gfx', 'pi_EITB1.jpg');
export_fig(fpath_fig)




%%
N = size(T, 1);
clc
clf
hold on
hs = [12];
Colors = {'b', 'r'};
LineStyles = {'-', '--', ':'};
LineWidths = [1, 1, 2];
for h = hs
    ma_UIAR = arrayfun(@(t) mean(UIAR(t:t+h)), 1:N-h)';
    ma_UITB = arrayfun(@(t) mean(UITB1(t:t+h)), 1:N-h)';
    Color = Colors{hs==h};
    LineStyle = LineStyles{hs==h};
    LineWidth = LineWidths(hs==h);
    plot(T.date(1+h:end), ma_UIAR,...
        'Color', Colors{1}, 'LineStyle', LineStyle, 'LineWidth', LineWidth)
    plot(T.date(1+h:end), ma_UITB,...
        'Color', Colors{2}, 'LineStyle', LineStyle, 'LineWidth', LineWidth)
    autocorr_(ma_UIAR, 1)
    autocorr_(ma_UITB, 1)
    std(ma_UIAR)
    std(ma_UIAR)/std(UIAR)

    std(ma_UITB)
    std(ma_UITB)/std(UITB1)
end

leg = {'UIAR', 'UITB'};
legend(leg, 'Location', 'north')

%title(sprintf('E(UITB)=%.3f', mean(UITB)))
%}

ylim([-.3 .3])

set(gcf, 'Color', 'w')
fpath_fig = fullfile(root, 'gfx', sprintf('UIAR_UITB_il_ma%d.jpg', hs));
export_fig(fpath_fig)





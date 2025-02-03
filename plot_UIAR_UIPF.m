% plot the ma and acf of UIAR and UIPF

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


min_date = datetime(2001, 4, 30); % comparing UIAR & UIPF
max_date = datetime(2019, 12, 31);
%max_date = datetime(2024, 3, 31);


li = 1;
fpath_UIPF = fullfile(lib_israel, sprintf('UIPF_%.1f.xlsx', li));
T = readtable(fpath_UIPF);
T.date = datetime(T.date);

idx = T.date >= min_date & T.date <= max_date;
T = T(idx, :);

UIPF = T.UIPF;

T.Pi = T.Pi_SA;
EIPF = T.Pi - UIPF;


p = 1;
d = 0;
q = 1;
forecast_lag = 2; % out-of-sample

fname_UIAR = sprintf('arima_%d%d%d_il.csv',p,d, q);
fpath_UIAR = fullfile(root, '..', 'NRC', fname_UIAR);

T = readtable(fpath_UIAR);
idx = T.date >= min_date & T.date <= max_date;
T = T(idx, :);


if forecast_lag > 0
    UIAR = T{:, sprintf('UIAR%d_%d%d%d', forecast_lag, p, d, q)};
else
    UIAR = T{:, sprintf('epsilon_%d%d%d', p, d, q)};
end

EIAR = T.Pi - UIAR;


N = size(T, 1);
hs = [36];
num_placebos = 1000;
placebos = std(UIAR)*randn(size(UIAR, 1), num_placebos);
placebo = placebos(:, 1);
ma_placebos = arrayfun(@(t) mean(placebos(t:t+hs, :), 1), 1:N-hs,...
    'UniformOutput', false)';
ma_placebos = cell2mat(ma_placebos);

p_std_wn = [...
    prctile(arrayfun(@(i) std(ma_placebos(:, i)), 1:size(ma_placebos, 2)), 2.5),...
    prctile(arrayfun(@(i) std(ma_placebos(:, i)), 1:size(ma_placebos, 2)), 97.5)];
p_wn = [...
    prctile(arrayfun(@(i) std(ma_placebos(:, i)), 1:size(ma_placebos, 2)), 2.5),...
    prctile(arrayfun(@(i) std(ma_placebos(:, i)), 1:size(ma_placebos, 2)), 97.5)];


idx = randi(size(UIAR, 1), [hs, 1000]);
UIAR_wr = UIAR(idx);
mean_UIAR_wr = mean(UIAR_wr, 1);
p_UIAR = [prctile(mean_UIAR_wr, 2.5), prctile(mean_UIAR_wr, 97.5)];



clf
subplot(1,2,1)
hold on
Colors = {'b', 'r', 'm', 'k'};
LineStyles = {'-', '-', '-', '--'};
LineWidths = [1, 1, 2];
for h = hs
    ma_UIAR = arrayfun(@(t) mean(UIAR(t:t+h)), 1:N-h)';
    ma_UIPF = arrayfun(@(t) mean(UIPF(t:t+h)), 1:N-h)';
    %ma_placebo = arrayfun(@(t) mean(placebo(t:t+h)), 1:N-h)';
    ma_placebo = ma_placebos(:, 1);
    Color = Colors{hs==h};
    LineStyle = LineStyles{hs==h};
    LineWidth = LineWidths(hs==h);
    plot(T.date(1+h:end), ma_UIAR,...
        'Color', Colors{1}, 'LineStyle', LineStyles{1}, 'LineWidth', LineWidth)

    
    plot(T.date(1+h:end), ma_UIPF,...
        'Color', Colors{3}, 'LineStyle', LineStyles{3}, 'LineWidth', LineWidth)
    %hline(p_wn, 'k--')

    %{
    plot(T.date(1+h:end), ma_placebo,...
        'Color', Colors{4}, 'LineStyle', LineStyles{4}, 'LineWidth', LineWidth)
    %}
    %autocorr_(ma_UIAR, 1)
    %autocorr_(ma_UIPF, 1)
    %autocorr_(placebo, 1)
    
    [std(UIAR), sqrt(mean(UIAR.^2)), std(ma_UIAR), std(ma_UIAR)/std(UIAR)]
    [std(UIPF), sqrt(mean(UIPF.^2)), std(ma_UIPF), std(ma_UIPF)/std(UIPF)]
    

end


leg = {'UIAR', 'UIPF', 'white noise'};
legend(leg, 'Location', 'north')


subplot(1,2,2)

histogram(std(ma_placebos))
vline(std(ma_UIAR), Colors{1})
vline(std(ma_UIPF), Colors{3})

corrcoef([UIAR, UIPF])


return
%set(gcf, 'Color', 'w')
%fpath_fig = fullfile(root, 'gfx', 'UIAR_UIPF_il_ma.jpg');
%export_fig(fpath_fig)

%%
clc; clf
%placebo = placebos(:, randi(size(placebos, 2)));
hs = 1:24;
acf_UIAR = arrayfun(@(h) autocorr_(UIAR, h), hs);
acf_UIPF = arrayfun(@(h) autocorr_(UIPF, h), hs);
acf_placebo = arrayfun(@(h) autocorr_(placebo, h), hs);
%
acf_placebos = arrayfun(@(i) arrayfun(@(h) autocorr_(placebos(:, i), h), hs),...
    1:size(placebos, 2), 'UniformOutput', false)';
acf_placebos = cell2mat(acf_placebos);
acf_placebo = median(acf_placebos, 1);
size(acf_placebos)


%}

%{
arrayfun(@(h) find(acf_UIAR(h) <= sort(acf_placebos(:, h)), 1, 'first')/...
    size(placebos, 2), hs)
arrayfun(@(h) find(acf_UIPF(h) <= sort(acf_placebos(:, h)), 1, 'first')/...
    size(placebos, 2), hs, 'UniformOutput', false)
%}

%[P,DW] = 


%Colors = {'b', 'r', 'm', 'k'}
bar(hs, [acf_UIAR; acf_UIPF; acf_placebo])

set(gcf, 'Color', 'w')
%fpath_fig = fullfile(root, 'gfx', 'UIAR_UIPF_il_acf.jpg');
%export_fig(fpath_fig)

leg = {'UIAR', 'UIPF', 'white noise'};
legend(leg)

%%


clf; clc
subplot(1,2,1)
hold on
plot(T.date, T.Pi, 'LineWidth', 1.5)
plot(T.date, EIAR, '--')
title('EIAR')


subplot(1,2,2)
hold on
plot(T.date, T.Pi, 'LineWidth', 1.5)
plot(T.date, EIPF, '--')
title('EIPF')


rmse(T.Pi, EIAR)
rmse(T.Pi, EIPF)

%%
clf; clc
subplot(2,1,1)
hold on
plot(T.date, EIAR)
plot(T.date, EIPF)
title('Expected inflation')
legend({'EIAR', 'EIPF'})
subplot(2,1,2)
hold on
plot(T.date, UIAR)
plot(T.date, UIPF)
title('Unexpected inflation')
legend({'UIAR', 'UIPF'})


[std(EIAR), std(EIPF)]
[std(UIAR), std(UIPF)]


%%

clf; clc
subplot(1,2,1)
hold on
plot(T.date, T.Pi, 'LineWidth', 1.5)
plot(T.date, EIAR, '--')
title('EIAR')


subplot(1,2,2)
hold on
plot(T.date, T.Pi, 'LineWidth', 1.5)
plot(T.date, EIPF, '--')
title('EIPF')


[rmse(T.Pi, EIAR), sqrt(mean(UIAR.^2))]
[rmse(T.Pi, EIPF), sqrt(mean(UIPF.^2))]

std(T.Pi)



%%

clf; clc
hold on
plot(T.date, EIAR, 'b-')
plot(T.date, EIPF, 'r-')
ylim([-1.5, +1.5])

legend({'EIAR', 'EIPF'})
set(gcf, 'Color', 'w')


[rmse(T.Pi, EIAR), sqrt(mean(UIAR.^2))]
[rmse(T.Pi, EIPF), sqrt(mean(UIPF.^2))]

std(T.Pi)


fpath_fig = fullfile(root, 'gfx', 'EIAR_EIPF_il.jpg');
export_fig(fpath_fig)



%%

clf; clc
hold on
plot(T.date, UIAR, 'b-')
plot(T.date, UIPF, 'r-')
ylim([-1.5, +1.5])

legend({'UIAR', 'UIPF', })
set(gcf, 'Color', 'w')


[rmse(T.Pi, EIAR), sqrt(mean(UIAR.^2))]
[rmse(T.Pi, EIPF), sqrt(mean(UIPF.^2))]

std(T.Pi)


fpath_fig = fullfile(root, 'gfx', 'UIAR_UIPF_il.jpg');
export_fig(fpath_fig)


%%

clf; clc
hold on
plot(T.date, EIAR, 'b-')
plot(T.date, EIPF, 'r-')
plot(T.date, T.Pi, 'k--')
ylim([-1.5, +1.5])

legend({'EIAR', 'EIPF', '\pi'})
set(gcf, 'Color', 'w')


[rmse(T.Pi, EIAR), sqrt(mean(UIAR.^2))]
[rmse(T.Pi, EIPF), sqrt(mean(UIPF.^2))]

std(T.Pi)


fpath_fig = fullfile(root, 'gfx', 'pi_EIAR_EIPF_il.jpg');
export_fig(fpath_fig)

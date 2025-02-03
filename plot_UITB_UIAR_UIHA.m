% plot the ma and acf of UITB, UIAR and UIHA

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


%min_date = datetime(1998, 2, 28); % comparing UITB & UIAR
min_date = datetime(2001, 10, 31); % comparing UITB & UIAR & UIHA
max_date = datetime(2024, 3, 31);


li = 1.5;
fpath_UIHA = fullfile(lib_israel, sprintf('UIHA_%.1f.xlsx', li));
T = readtable(fpath_UIHA);
T.date = datetime(T.date);

idx = T.date >= min_date & T.date <= max_date;
T = T(idx, :);

UIHA = T.UIHA;

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



%fname_UITB = 'MAKAM_yields_M01.M.xlsx';
fname_UITB = 'BL.TELBOR_01M.M.xlsx';
fpath_UITB = fullfile(root, 'ML_FGLS_il_ZRD_oos', 'Rf', fname_UITB);

T = readtable(fpath_UITB);
T.date = T.date + 1;
T.date = T.date - calmonths(1);
T.date = T.date - 1;


idx = T.date >= min_date & T.date <= max_date;
T = T(idx, :);

%UITB = T.err_t_1__oos; % beta estimated
UITB = T.err_1_t_1__oos; % beta(1) = 1
Pi = T.Pi_t_1_;

N = size(T, 1);
hs = [36];
num_placebos = 1000;
placebos = std(UIAR)*randn(size(UIAR, 1), num_placebos);
placebo = placebos(:, 1);
ma_placebos = arrayfun(@(t) mean(placebos(t:t+hs, :), 1), 1:N-hs,...
    'UniformOutput', false)';
ma_placebos = cell2mat(ma_placebos);

p_wn = [...
    prctile(arrayfun(@(i) std(ma_placebos(:, i)), 1:size(ma_placebos, 2)), 2.5),...
    prctile(arrayfun(@(i) std(ma_placebos(:, i)), 1:size(ma_placebos, 2)), 97.5)];


idx = randi(size(UITB, 1), [hs, 1000]);
UITB_wr = UITB(idx);
mean_UITB_wr = mean(UITB_wr, 1);
p_UITB = [prctile(mean_UITB_wr, 2.5), prctile(mean_UITB_wr, 97.5)];


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
    ma_UITB = arrayfun(@(t) mean(UITB(t:t+h)), 1:N-h)';
    ma_UIHA = arrayfun(@(t) mean(UIHA(t:t+h)), 1:N-h)';
    %ma_placebo = arrayfun(@(t) mean(placebo(t:t+h)), 1:N-h)';
    ma_placebo = ma_placebos(:, 1);
    Color = Colors{hs==h};
    LineStyle = LineStyles{hs==h};
    LineWidth = LineWidths(hs==h);
    plot(T.date(1+h:end), ma_UIAR,...
        'Color', Colors{1}, 'LineStyle', LineStyles{1}, 'LineWidth', LineWidth)
    plot(T.date(1+h:end), ma_UITB,...
        'Color', Colors{2}, 'LineStyle', LineStyles{2}, 'LineWidth', LineWidth)
    
    hline(p_wn, 'k--')
    %{
    plot(T.date(1+h:end), ma_UIHA,...
        'Color', Colors{3}, 'LineStyle', LineStyles{3}, 'LineWidth', LineWidth)
    %}
    %{
    plot(T.date(1+h:end), ma_placebo,...
        'Color', Colors{4}, 'LineStyle', LineStyles{4}, 'LineWidth', LineWidth)
    %}
    %autocorr_(ma_UIAR, 1)
    %autocorr_(ma_UITB, 1)
    %autocorr_(ma_UIHA, 1)
    %autocorr_(placebo, 1)
    
    [std(UIAR), sqrt(mean(UIAR.^2)), std(ma_UIAR), std(ma_UIAR)/std(UIAR)]
    [std(UITB), sqrt(mean(UITB.^2)), std(ma_UITB), std(ma_UITB)/std(UITB)]
    [std(UIHA), sqrt(mean(UIHA.^2)), std(ma_UIHA), std(ma_UIHA)/std(UIHA)]
    

end

corrcoef([UIAR, UITB, UIHA])


leg = {'UIAR', 'UITB', 'UIHA', 'white noise'};
legend(leg, 'Location', 'north')

%title(sprintf('E(UITB)=%.3f', mean(UITB)))
%}

subplot(1,2,2)

histogram(std(ma_placebos))
vline(std(ma_UIAR), Colors{1})
vline(std(ma_UITB), Colors{2})
vline(std(ma_UIHA), Colors{3})

return
%set(gcf, 'Color', 'w')
%fpath_fig = fullfile(root, 'gfx', 'UIAR_UITB_UIHA_il_ma.jpg');
%export_fig(fpath_fig)

%%
clc; clf
%placebo = placebos(:, randi(size(placebos, 2)));
hs = 1:24;
acf_UIAR = arrayfun(@(h) autocorr_(UIAR, h), hs);
acf_UITB = arrayfun(@(h) autocorr_(UITB, h), hs);
acf_UIHA = arrayfun(@(h) autocorr_(UIHA, h), hs);
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
arrayfun(@(h) find(acf_UITB(h) <= sort(acf_placebos(:, h)), 1, 'first')/...
    size(placebos, 2), hs)
arrayfun(@(h) find(acf_UIHA(h) <= sort(acf_placebos(:, h)), 1, 'first')/...
    size(placebos, 2), hs, 'UniformOutput', false)
%}

%[P,DW] = 


%Colors = {'b', 'r', 'm', 'k'}
bar(hs, [acf_UIAR; acf_UITB; acf_UIHA; acf_placebo])

set(gcf, 'Color', 'w')
%fpath_fig = fullfile(root, 'gfx', 'UITB_il_acf.jpg');
%export_fig(fpath_fig)

leg = {'UIAR', 'UITB', 'UIHA', 'white noise'};
legend(leg)

%%

EIAR = Pi - UIAR;
EITB = Pi - UITB;
EIHA = Pi - UIHA;

clf; clc
subplot(2,2,1)
hold on
plot(T.date, Pi, 'LineWidth', 1.5)
plot(T.date, EIAR, '--')

subplot(2,2,2)
hold on
plot(T.date, Pi, 'LineWidth', 1.5)
%plot(T.date(1:end-1), EITB(2:end))
plot(T.date, EITB, '--')

rmse(Pi, EIAR)
rmse(Pi, EITB)
rmse(Pi, EIHA)


subplot(2,2,3)
hold on
plot(T.date, Pi, 'LineWidth', 1.5)
plot(T.date, EIHA, '--')








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


fname_UITB = 'MAKAM_yields_M01.M.xlsx';
fpath_UITB = fullfile(root, 'ML_FGLS_il_ZRD_oos', fname_UITB);

T = readtable(fpath_UITB);
T.date = T.date + 1;
T.date = T.date - calmonths(1);
T.date = T.date - 1;

min_date = datetime(1998, 2, 28);
max_date = datetime(2024, 3, 31);
idx = T.date >= min_date & T.date <= max_date;
T = T(idx, :);

UITB = T.err_t_1__oos;

N = size(T, 1);

clf
hold on
hs = [12, 24, 36];
Colors = {'k', 'b', 'r'};
LineStyles = {'-', '--', ':'};
LineWidths = [1, 1, 2];
for h = hs
    ma = arrayfun(@(t) mean(UITB(t:t+h)), 1:N-h)';
    Color = Colors{hs==h};
    LineStyle = LineStyles{hs==h};
    LineWidth = LineWidths(hs==h);
    plot(T.date(1+h:end), ma,...
        'Color', Color, 'LineStyle', LineStyle, 'LineWidth', LineWidth)
    autocorr_(ma, 1)
end

leg = arrayfun(@num2str, hs, 'UniformOutput', false);
leg = strcat({'ma-'}, leg, {' months'});
legend(leg, 'Location', 'north')

%title(sprintf('E(UITB)=%.3f', mean(UITB)))
%}

set(gcf, 'Color', 'w')
fpath_fig = fullfile(root, 'gfx', 'UITB_il_ma.jpg');
export_fig(fpath_fig)

return
%%
clc; clf
hs = 1:12;
acf_UITB = arrayfun(@(h) autocorr_(UITB, h), hs);

arrayfun(@(h) autocorr_(UITB, h), hs)
bar(hs, [acf_UITB; ])

set(gcf, 'Color', 'w')
fpath_fig = fullfile(root, 'gfx', 'UITB_il_acf.jpg');
export_fig(fpath_fig)











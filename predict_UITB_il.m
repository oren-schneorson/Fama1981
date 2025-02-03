
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



%series_Rf = 'MAKAM_yields_M01.M_15';
series_Rf = 'MAKAM_yields_M01.M';
%series_Rf = 'BL.TELBOR_01M.M';
if strcmp(series_Rf, 'BL.TELBOR_01M.M')
    min_date = datetime(2001, 10, 31);
    max_date = datetime(2019, 12, 31);
elseif strcmp(series_Rf, 'MAKAM_yields_M01.M')
    min_date = datetime(2000, 1, 31);
    max_date = datetime(2019, 12, 31);
else
    min_date = datetime(2000, 1, 31);
    max_date = datetime(2019, 12, 31);
end


fname_UITB = [series_Rf, '.xlsx'];
fpath_UITB = fullfile(root, 'ML_FGLS_il_ZRD_oos', 'Rf', fname_UITB);

T = readtable(fpath_UITB);
%{
T.date = T.date + 1;
T.date = T.date - calmonths(1);
T.date = T.date - 1;
%}
idx = T.date >= min_date & T.date <= max_date;
T = T(idx, :);


fpath_fear_index = fullfile(lib_israel, 'FEAR_INDEX', 'M', 'FEAR_INDEX.csv');
FEAR_INDEX = readtable(fpath_fear_index);

T = innerjoin(T, FEAR_INDEX, 'Keys', 'date');


UITB = T.err_t_1__oos;
UITB1 = T.err_1_t_1__oos;
Pi_TB = T.Pi_t_1_;

Y = UITB1;
X = T.FEAR_INDEX;
X = [ones(size(X, 1), 1), X];

[B, BINT, R, ~, STATS] = regress(Y, X);
[BINT(:, 1), B, BINT(:, 2)]
STATS(1)

hs = 1:36;
arrayfun(@(h) autocorr_(R, h), hs)
arrayfun(@(h) autocorr_(Y, h), hs)



















%{
Fear index, Yossi Saadon
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

fpath = fullfile(lib_israel, 'FEAR_INDEX', 'D', 'FEAR_INDEX.csv');
FEAR_INDEX = readtable(fpath);



%%

clf; clc
plot(FEAR_INDEX.date, FEAR_INDEX.FEAR_INDEX)
p = prctile(FEAR_INDEX.FEAR_INDEX, 95);
hline(p)









%{


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


Range = 'A11:B5050'; % read enough rows to cover > (2022-1982) * 12
Range_meta_data = 'B8:B8';

datevar = 'observation_date';

datasets = {...
    'GDP_real',...
    'GDP',...
    'GDP_constant_per_capita',...
    'GDP_real_constant_national_prices',...
    'GDP_growth',...
    };
datasets = {...
    %'expected_inflation_TIPS',...
    'expected_inflation_Cleveland_model',...
    };

for dataset = datasets

    dataset = char(dataset);
    lib_series = fullfile(lib_data, dataset);

    d = dir(lib_series);
    d = d(~[d.isdir]);

    counter = 0;
    meta_data = table;


    for fname = {d.name}

        counter = counter + 1;
        fname = char(fname);
        series = strsplit(fname, '.');
        series = series(1);

        path = fullfile(lib_series, fname);
        meta_data_ = readtable(path, 'Range', Range_meta_data,...
            'ReadVariableNames', false);
        meta_data_ = meta_data_{1,1};
        meta_data_ = strsplit(char(meta_data_), ', ');

        desc = meta_data_(2);
        country = {'United States'};
        unit = meta_data_(2);
        freq = meta_data_(3);
        seasonal_adj = meta_data_(4);

        meta_data_ = struct;
        meta_data_.desc = desc;
        meta_data_.country = strrep(country, 'the ', '');
        meta_data_.unit = unit;
        meta_data_.freq = freq;
        meta_data_.seasonal_adj = ~contains(seasonal_adj, 'Not');
        meta_data_.series = series;

        T_ = readtable(path, 'Range', Range);

        idx = ~isnan(T_.(char(series)));

        if idx(end) == 1
            error('Range not enough')
        end

        T_ = T_(idx, :);

        meta_data_.start_date = datestr(T_{1, datevar}, 'yyyy-mm-dd');
        meta_data_.end_date = datestr(T_{end, datevar}, 'yyyy-mm-dd');

        meta_data = [meta_data; struct2table(meta_data_)];


        if ~exist('T', 'var')    
            T = T_;
        else
            T = outerjoin(T, T_, 'Keys', {datevar},...
                'MergeKeys', true);
        end


        fprintf('%s -- %.2f%%\n', char(country), counter/numel(d)*100)


    end


    %{
    T.(datevar) = cellfun(@(c) [c(7:10), c(4:5), c(1:2)],...
        T.(datevar), 'UniformOutput', false);
    %}
    T = sortrows(T, datevar);
    T.(datevar) = datestr(T.(datevar), 'yyyy-mm-dd');


    %{
    T.(datevar) = cellfun(@(c) [c(1:4), '-', c(5:6), '-', c(7:8)],...
        T.(datevar), 'UniformOutput', false);
    %}



    pat = '\d{1,2}';
    vars = T.Properties.VariableNames(2:end);
    y = cellfun(@(c) regexp(c,pat,'match'), vars);
    y = str2double(y);
    [~, col_idx] = sort(y);
    meta_data = meta_data(col_idx,:);
    col_idx = [1, 1+col_idx];
    T = T(:,col_idx);

    path = fullfile(lib_data, [dataset, '.xlsx']);
    writetable(T, path)

    path = fullfile(lib_data, [dataset, '_meta_data.xlsx']);
    writetable(meta_data, path)
    clear T
    clear meta_data


end

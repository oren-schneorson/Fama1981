%{

This script replicates desc stats in tables 1-2, Nelson (1976)+
extension to include breakdown to expected and unexpected inflation.

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
matlab_dir_TASE = fullfile(matlab_dir, 'TASE');

lib_data = fullfile('/media', username, 'D', 'data');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')

lib_cpi = fullfile(lib_data, 'cpi_usa');
lib_EITB = fullfile(root, 'EITB');
lib_UITB = fullfile(root, 'UITB');
lib_crsp = fullfile(lib_data, 'fama_french_factors', 'idx_mth');
lib_ir = fullfile(lib_data, 'ir');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');



% flag to override period in downstream scripts: false by default
lag_inflation = 1;
log_approx = true;
%log_approx = false;
freq = 'M';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
to_Quarterly = to_Annual;
vintage = true;

pcent_mult = 1; % set as 100 for percent

%{
% original dates: Nelson 1976
sample_tab = table(...
    [1:4]',...
    [datetime(1953,6,30), datetime(1953,6,30), datetime(1953,6,30), datetime(1964,1,31)]',...
    [datetime(1971,4,30), datetime(1974,2,28), datetime(1963,12,31), datetime(1974,2,28)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%{
% extended dates
sample_tab = table(...
    [1:8]',...
    [datetime(1953,1,31), datetime(1953,1,31), datetime(2000,1,31),...
     datetime(1953,1,31), datetime(2007,7,31), datetime(1953,12,31),...
     datetime(1954,1,31), datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(2019,12,31), datetime(1977,12,31),...
    datetime(1976,12,31), datetime(1977,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%
% paper dates
sample_tab = table(...
    [1:4]',...
    [datetime(1953,1,31), datetime(1953,1,31),...
     datetime(1953,1,31), datetime(2000,1,31)]',...
    [datetime(2019,12,31), datetime(2022,12,31),...
     datetime(1999,12,31), datetime(2019,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

tab_EI = table;
tab_UI = table;
tab_EI_UI = table;

VariableNames_EI = {'period', 'min_date', 'max_date',...
    'constant', 'beta_EI',...
    'tstat_constant', 'tstat_beta_EI',...
    'DW', 'R_sq', 'R_sq-adj', 'F', 'F_nc'};
VariableNames_UI = {'period', 'min_date', 'max_date',...
    'constant', 'beta_UI',...
    'tstat_constant', 'tstat_beta_UI',...
    'DW', 'R_sq', 'R_sq-adj', 'F', 'F_nc'};
VariableNames_EI_UI = {'period', 'min_date', 'max_date',...
    'constant', 'beta_EI', 'beta_UI',...
    'tstat_constant', 'tstat_beta_EI', 'tstat_beta_UI',...
    'DW', 'R_sq', 'R_sq-adj', 'F', 'F_nc'};
RowNames = {'Const.', 'beta_EI', 'beta_UI'};


series_CPI = 'CPIAUCSL';
%series_CPI = 'CWUR0000SA0';

series_Rf = 'RF_FF';
series_crsp = 'CRSP_vw';





data = struct(...
    'series', {...
    'CWUR0000SA0',... CPI, not seasonally adjusted, BLS
    'CPIAUCSL',... Price index of all urban consumers
    series_crsp,...
    series_Rf,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_crsp,...
    lib_ir,...
    %'',...
    },...
    'src', {...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    %'',...
    },...
    'VariableNames', cell(1,4)...
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
        %{
        if strcmp(data(i).series, 'CWUR0000SA0')
            cpi_nsa = T_.CWUR0000SA0;
            cpi_nsa(T_.date.Year< 1947) = NaN;
            cpi_nsa(T_.date.Year> 1979) = NaN;
            cpi_sa = sa_adj(cpi_nsa, 12);
            T_.CWUR0000SA0 = cpi_sa;
        end
        %}
        %{
        clf
        hold on
        plot(T_.date, log(cpi_nsa))
        plot(T_.date, log(cpi_sa))
        return
        %}
        
    end

    
    if ~exist('T', 'var')
        T = T_;
        metadata = metadata_;
    else
        if any(ismember(T.Properties.VariableNames, data(i).series))
            continue
        end

        %{
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true, 'Type', 'left');
        %}
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true');
        metadata = [metadata; metadata_];
        
    end
    

end


fpath = fullfile(matlab_dir, 'NRC', 'arima_101_usa.csv')

arima_101 = readtable(fpath);
arima_101 = renamevars(arima_101, 'Pi_hat', 'Pi_hat_101');

fpath = fullfile(root, 'TB_model_usa.csv')
TB_model = readtable(fpath);
TB_model = renamevars(TB_model, 'max_date', 'date');
TB_model = renamevars(TB_model, 'Pi_hat', 'Pi_hat_TB');
TB_model = renamevars(TB_model, 'err', 'epsilon_TB');

if strcmp(freq, 'M')
    %pass
elseif strcmp(freq, 'A')
    T = table2timetable(T);
    T = retime(T, 'yearly', to_Annual);
    T = timetable2table(T);
elseif strcmp(freq, 'Q')
    T = table2timetable(T);
    T = retime(T, 'quarterly', to_Quarterly);
    T = timetable2table(T);
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

%T = fillmissing(T, 'previous');

T = outerjoin(T, arima_101(:, {'date', 'Pi', 'Pi_hat_101', 'epsilon_101'}),...
    'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
T = outerjoin(T, TB_model(:, {'date', 'Pi_hat_TB', 'epsilon_TB'}),...
    'type', 'Left', 'Keys', 'date', 'MergeKeys', true);




% Sample period
%for period = sample_tab.period(end-2)'
for period = sample_tab.period(2)'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);
    fprintf('%s-%s\n', min_date, max_date)
    
    %{
    series_EITB = sprintf('EITB__%s-%s.M',...
        min_date, max_date);
    series_UITB = sprintf('UITB__%s-%s.M',...
        min_date, max_date);
    
    fname = [series_EITB, '.csv'];
    fpath_EITB = fullfile(lib_EITB, fname);
    EITB = readtable(fpath_EITB);

    fname = [series_UITB, '.csv'];
    fpath_UITB = fullfile(lib_UITB, fname);
    UITB = readtable(fpath_UITB);

    T = outerjoin(T, EITB,...
        'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
    T = outerjoin(T, UITB,...
        'type', 'Left', 'Keys', 'date', 'MergeKeys', true);
    
    %}
idx = ~isnan(T.epsilon_TB) & ~isnan(T.epsilon_101);
clf
hold on
%plot(T.date(idx), T.Pi_hat_TB(idx)+T.epsilon_TB(idx))
%plot(T.date(idx), T.Pi_hat_101(idx)+T.epsilon_101(idx), '--r')
%plot(T.date(idx), T.Pi(idx), 'k-')

LineWidth = 1;
plot(T.date(idx), T.Pi(idx), 'r:')
plot(T.date(idx), T.Pi_hat_TB(idx), 'k-', 'LineWidth', LineWidth)
plot(T.date(idx), T.Pi_hat_101(idx), 'k--', 'LineWidth', LineWidth)

%plot(T.date(idx), T.epsilon_TB(idx), 'b-')
%plot(T.date(idx), T.epsilon_101(idx), '-r')

%plot(T.date(idx), abs(T.epsilon_TB(idx))-abs(T.epsilon_101(idx)))
%hline(mean(abs(T.epsilon_TB(idx))-abs(T.epsilon_101(idx))))

%{
T.Pi(abs(T.Pi)<1e-3) = NaN;

plot(T.date(idx), T.epsilon_TB(idx)./T.Pi(idx))
plot(T.date(idx), T.epsilon_101(idx)./T.Pi(idx))
%}

%plot(T.date(idx), T.EITB(idx))
%plot(T.date(idx), T.Pi_hat_TB(idx))

%plot(T.date(idx), T.epsilon_101(idx))
%plot(T.date(idx), T.epsilon_TB(idx))

xlabel('Date')
ylabel('Percent')

legend({'CPI inflation', 'ARMA(1,1)', 'TBill'})
%legend({'ARMA(1,1)', 'TBill'})

set(gcf, 'Position', [169 225 931 595]) % [width=0.88\textwidth]
set(gca, 'FontSize', 16)
set(gca, 'FontName', 'Lucida')

set(gcf, 'Color', 'w')

fpath_fig = fullfile(root, 'gfx', 'comparing_inflation_prediction_models_usa.png');
export_fig(gcf, fpath_fig)

min(T.date(idx))
max(T.date(idx))
sqrt(mean(T.epsilon_101(idx).^2, 'omitnan'))
sqrt(mean(T.epsilon_TB(idx).^2, 'omitnan'))
%sqrt(mean(T.UITB(idx).^2, 'omitnan'))
std(T.Pi(idx))

corr(T.epsilon_101(idx), T.epsilon_TB(idx))
corr(T.Pi_hat_101(idx), T.Pi_hat_TB(idx))





















    
end


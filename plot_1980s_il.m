%{

This script runs replicates desc stats in table 4, Fama (1981).

It does so with updated data, which is a bit different to the original,
and also with vintage data from the Survey of Current Business (SCB) of the
Bureau of Economic Analysis (BEA).

For example, some series have had revisions. I did the best I could to look
at the original data, even going to the Survey of Current Business real
time publications to see some of the series Fama used. In the end I decided
against using those series, because results are not very different when
using the updated series.

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

lib_data = fullfile('/media', username, 'D', 'data', 'BOI');
lib_israel = fullfile(lib_data, 'Israel');

addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')


% flag to override period in downstream scripts: false by default
lag_inflation = 0;
log_approx = true;
%log_approx = false;
freq = 'M';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
vintage = true;

pcent_mult = 1; % set as 100 for percent

sample_tab = table(...
    [1:6]',...
    [...
    datetime(1980,1,31), datetime(1980,1,31), datetime(1986,1,31),...
    datetime(1980,1,31), datetime(1999,12,31), datetime(2008,10,31),...
     ]',...
    [...
    datetime(2024,6,31), datetime(2008,9,30), datetime(2008,9,30),...
    datetime(1985,12,31), datetime(2008,9,30), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});


table4 = table;

series_CPI = 'CP_SA.M';
series_CPI_NSA = 'CP_NSA.M';
EI_measure = {'AR'};

if strcmp(freq, 'M')
    series_Rf = 'TSB_BAGR_MAKAM_01M.M'; % 1 month MAKAM, annualized
elseif strcmp(freq, 'Q')
    series_Rf = 'TSB_BAGR_MAKAM_03M.M'; % 3 month MAKAM, annualized
elseif strcmp(freq, 'A')
    series_Rf = 'TSB_BAGR_MAKAM_12M.M'; % 12 month MAKAM, annualized
else
    error('freq must be M, Q or A')
end

series_ms = 'ZA108.M'; % waiting for Elad from MOS
%series_ms = {'A346.M_E', 'ZA215.M'}; % ask Elad about RESNALNS, 

series_ea = 'mdd13.level.m'; % excl diamonds
%series_ea = 'CLS11.TPR_C.Q_SA_CHAINED'; % excl diamonds


%series_cx = 'AM_INV_TOTALT_Q_N';
series_cx = 'AM_INV_ISKI_Q_N';

%series_gnp = 'GNP.Q_N';
series_gnp = 'GNP.Q_FP';
% residential and nonresidential equipment and structures, capital stock
series_ns = {'AM_NSTOCK_ISKI_Q_N'}; % includes inventories... not exactly like Fama
series_gs = {'AM_NSTOCK_ISKI_Q_N'}; % includes inventories... not exactly like Fama
series_roc = 'ROC_BS.Q';
series_roc_denom = {'AM_NSTOCK_TOTALT_Q_N'}; % series_ns_fred + inventories
series_tase = 'SPASTT01ILM661N';
series_C = {'C_@DUR.Q_FP_SA'}; % excl durable, should contain services


lib_cpi = fullfile(lib_israel);
lib_ms = fullfile(lib_israel, 'monagg');
lib_mbase = fullfile(lib_israel, 'money_base');
lib_reserves = fullfile(lib_israel, 'reserves');
lib_roc = fullfile(lib_israel);


lib_gnp = fullfile(lib_israel, 'GNP');  % gross national product
lib_capital_stock = fullfile(lib_israel, 'capital_stock');
lib_cx = fullfile(lib_capital_stock, 'inv');  % capital expenditure (investment)
lib_ns = fullfile(lib_capital_stock, 'nstock');  % net capital stock
lib_gs = fullfile(lib_capital_stock, 'gstock');  % gross capital stock

%lib_ea = fullfile(lib_israel, 'indprod');  % economic activity
lib_ea = fullfile(lib_israel);  % economic activity, monthly series

%lib_ir = fullfile(lib_israel);
lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M', 'with_metadata');
lib_pop = fullfile(lib_israel, 'POP');
lib_consumption = fullfile(lib_israel, 'C');
%lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');

lib_stock_markets = fullfile(lib_data, 'stock_markets');

lib_EITB = fullfile(root, 'EITB_il');
lib_UITB = fullfile(root, 'UITB_il');



data = struct(...
    'series', {...
    series_CPI,... General Consumer Price Index, seasonally adjusted
    series_CPI_NSA,... General Consumer Price Index, non seasonally adjusted
    series_Rf,...
    'ZA108.M',...
    'A346.M_E',...
    'ZA215.M',...
    series_gnp,... % real GNP
    series_ea,...
    series_cx,...
    'AM_NSTOCK_ISKIM_Q_N',... % Business Sector; Equipment
    'AM_NSTOCK_ISKIS_Q_N',... % Business Sector; Structures
    'AM_NSTOCK_BILDRES_Q_N',...  % Construction Sector; Residential All
    'AM_NSTOCK_TOTALT_Q_N',...
    'AM_NSTOCK_ISKI_Q_N',...
    'AM_GSTOCK_ISKIM_Q_N',... % Business Sector; Equipment
    'AM_GSTOCK_ISKIS_Q_N',... % Business Sector; Structures
    'AM_GSTOCK_BILDRES_Q_N',...  % Construction Sector; Residential All
    'AM_GSTOCK_TOTALT_Q_N',...
    'AM_GSTOCK_ISKI_Q_N',...
    series_roc,...
    'SPASTT01ILM661N',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_ir,...
    lib_mbase,...
    lib_mbase,...
    lib_reserves,...
    lib_gnp,...
    lib_ea,...
    lib_cx,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_ns,...
    lib_gs,...
    lib_gs,...
    lib_gs,...
    lib_gs,...
    lib_gs,...
    lib_roc,...
    lib_stock_markets,...
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
    'FAME',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'NONE',...
    'FAME',...
    'FRED',...
    },...
    'VariableNames', cell(1,21)...
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
        %metadata = [metadata; metadata_(:, metadata.Properties.VariableNames)];
        
    end
    
end

for period = 1:4

min_date = sample_tab.min_date(period);
max_date = sample_tab.max_date(period);

idx_sample = true(size(T, 1), 1);
idx_sample = idx_sample & T.date >= min_date;
idx_sample = idx_sample & T.date <= max_date;    

idx_col_CPI_NSA = ismember({data.series}, series_CPI_NSA);


x = T.(series_tase)(idx_sample);
y = T{idx_sample, [data(idx_col_CPI_NSA).VariableNames]};
z = x./y;

x = log(x/x(1));
y = log(y/y(1));
z = log(z/z(1));



clf
set(gcf, 'Color', 'w')
set(gca, 'Fontsize', 14)
set(gcf, 'Position', [71   189   915   683])

hold on
plot(T.date(idx_sample), x, 'b-')
plot(T.date(idx_sample), y, 'r-')

yyaxis right
plot(T.date(idx_sample), z, 'k--')
ax = gca;
ax.YAxis(2).Color = 'k';
ylim([0, ax.YLim(2)])

legend({'Log(TASE-all)', 'Log(CPI) (non seasonally adjusted)', 'Real TASE index (right)'},...
    'Location', 'northwest')

fname_fig = sprintf('stocks_inflation_il__%s--%s.png', char(min_date), char(max_date));
fpath_fig = fullfile(root, 'gfx', fname_fig);
export_fig(fpath_fig)

end


    



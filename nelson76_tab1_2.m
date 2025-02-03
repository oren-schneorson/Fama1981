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
lag_inflation = 0;
pcent_mult = 1; % set as 100 for percent
log_approx = true;
%log_approx = false;
freq = 'M';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
to_Quarterly = to_Annual;
vintage = true;

%{
if ~vintage
    lib_ns_vintage = lib_ns;
end
%}


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

VariableNames1 = {'period', 'min_date', 'max_date',...
    'lead_lag', 'B', 't_stat', 'corr'};
VariableNames2 = {'period', 'min_date', 'max_date',...
    'constant', 'lag0', 'lag1', 'lag2', 'lag3', 'lag4','sum_lags',...
    't_stat_constant', 't_stat_lag0', 't_stat_lag1', 't_stat_lag2', 't_stat_lag3', 't_stat_lag4','t_stat_sum_lags',...
    'DW', 'R_sq', 'F', 'F_nc'};


tab1 = table;
tab2 = table;
RowNames2 = arrayfun(@(lag) sprintf('lag%d', lag),...
    0:4, 'UniformOutput', false);


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

% Sample period
%for period = sample_tab.period(end-2)'
for period = sample_tab.period'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);
    
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

    T = outerjoin(T, EITB, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    T = outerjoin(T, UITB, 'type', 'Left', 'Keys', 'date',...
        'MergeKeys', true);
    
    EITB = T.EITB;
    UITB = T.UITB;
    T = removevars(T, {'EITB', 'UITB'});  

    if strcmp(freq, 'A')
        EITB = EITB*12;
        UITB = UITB*12;
    elseif strcmp(freq, 'Q')
        EITB = EITB*4;
        UITB = UITB*4;
    
    elseif strcmp(freq, 'M')
        % pass
    else
        error('freq must be A, Q or M')
    end
        
    % ***************************************
    % CRSP value weighted inedx

    idx = ismember({data.series}, series_crsp);
    CRSP = T{:, data(idx).VariableNames};
    
    % ***************************************
    % Risk free rate

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames};

    % ***************************************
    % CPI

    idx = ismember({data.series}, series_CPI);
    CPI = T{:, data(idx).VariableNames};
    Pi = [NaN;...
        T{2:end,   data(idx).VariableNames}./...
        T{1:end-1, data(idx).VariableNames}];
    
    
    dCRSP = [NaN;...
        CRSP(2:end, :)./...
        CRSP(1:end-1, :)];
    
   
    if log_approx
        
        %Pi = log(Pi)*pcent_mult;
        dCRSP = log(dCRSP)*pcent_mult;

    else

        %Pi = (Pi - 1)*pcent_mult;
        dCRSP = (dCRSP- 1)*pcent_mult;

    end
    
    Pi = (Pi - 1)*pcent_mult;
    RS = dCRSP-Pi; % real stock return
    
    EITB = EITB/100*pcent_mult;
    UITB = UITB/100*pcent_mult;
    
    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;    

    Y = dCRSP;
    X2 = Pi;
    
    % regressions for Table 2
    X2 = arrayfun(@(lag) [...
        NaN(lag, size(X2,2)); X2(1:end-lag,:)...
        ], 0:1:4, 'UniformOutput', false);
    X2 = cell2mat(X2);
    K = size(X2, 2);

    RowNames2 = [{'Constant'}, RowNames2];
    X2 = [ones(size(X2, 1), 1), X2];

    Y = Y(idx_sample, :);
    X2 = X2(idx_sample, :);
    
    N = size(X2, 1);
    K = size(X2, 2);

    [B,BINT,R,~,STATS] = regress(Y, X2);
    
    sigma_hat = sqrt(R'*R/(N-K));

    H = [0, ones(1,5)];
    H_nc = eye(K);
    H_nc = H_nc(2:end,:); % for F without constant
    D = X2'*X2; % information matrix
    
    size(H_nc)
    size(inv(X2'*X2))
    
    D_nc = inv(H_nc * inv(D) * H_nc');

    F = 1/(sigma_hat^2 * K)*...
        B'*D*B;

    %{
        % matlab computes F_nc, not F
    STATS(2)
    F_nc = 1/(sigma_hat^2 * size(H_nc, 1))*...
        (H_nc*B)'*D_nc*(H_nc*B)
        %}


    V_hat = sigma_hat^2*inv(D);
    STDB = sqrt(diag(V_hat));
    t_stat_sum_B = H*B/sqrt(H*V_hat*H');


    DW = sum(diff(R).^2)/sum(R.^2); % Durbin-Watson stat


    t_stat_B = B./STDB;
    t_stat_B = mat2cell(t_stat_B, ones(size(t_stat_B)), size(t_stat_B,2));

    B_cell = mat2cell(B, ones(size(B)), size(B,2));

    tab_row = table(...
        period, min(T.date(idx_sample)), max(T.date(idx_sample)),...
        B_cell{:}, H*B, t_stat_B{:}, t_stat_sum_B, DW, STATS(1), F, STATS(2),...
        'VariableNames', VariableNames2);
    tab2 = [tab2; tab_row];
    
    
    for lead_lag = -4:4

    Y = dCRSP;
    X1 = Pi;

    if lead_lag < 0
        X1 = [NaN(abs(lead_lag), size(X1,2)); X1(1:end-abs(lead_lag),:)];
    else
        X1 = [X1(max(1, lead_lag+1):end,:); NaN(lead_lag, size(X1,2))];
    end    

    RowNames1 = {'Constant', 'Inflation'};
    X1 = [ones(size(X1, 1), 1), X1];


    fprintf('Full sample, %s--%s; lag=%d, lead=%d\n',...
        min(T.date(idx_sample)), max(T.date(idx_sample)), abs(lead_lag), lead_lag)
    fprintf('PI between %.3f-%.3f, Aveage inflation: %.3f\n',...
        prctile(Pi(idx_sample), 0),...
        prctile(Pi(idx_sample), 100),...
        mean(X1(:,end), 'omitnan'))

    X1 = X1(idx_sample, :);
    Y = Y(idx_sample);

    N = size(X1, 1);
    K = size(X1, 2);
    
    
    % Full
    [B,BINT,~,~,STATS] = regress(Y, X1);
    
    table(BINT(:,1), B, BINT(:,2), 'VariableNames',...
        {'Lower', 'Coef', 'Upper'}, 'RowNames',...
        RowNames1)
    [B_, STDB_] = lscov(X1, Y);

    tab_row = table(...
        period, min(T.date(idx_sample)), max(T.date(idx_sample)), lead_lag,...
        B(end), B_(end)./STDB_(end), corr(X1(:, end), Y),...
        'VariableNames', VariableNames1);
    tab1 = [tab1; tab_row];

    end

    
end


PycharmProjects_dir = fullfile('/home', username, 'PycharmProjects');

fname = 'data_tab1.csv';
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', fname);
writetable(tab1, fpath)

fname = 'data_tab2.csv';
fpath = fullfile(PycharmProjects_dir, 'Stocks_inflation', 'data', fname);
writetable(tab2, fpath)


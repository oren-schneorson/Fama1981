%{

This script replicates desc stats in tables 1-2, Nelson (1976).

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
lib_data = fullfile('/media', username, 'D', 'data', 'BOI');

lib_israel = fullfile(lib_data, 'Israel');
lib_cpi = lib_israel;


addpath(fullfile(matlab_dir, 'altmany-export_fig-410f0ad')
addpath(fullfile(matlab_dir, 'my_functions')


% flag to override period in downstream scripts: false by default
lag_inflation = 0;
log_approx = true;
%log_approx = false;
freq = 'M';
to_Annual = 'lastvalue';
%to_Annual = 'mean';
to_Quarterly = to_Annual;
vintage = true;
pcent_mult = 100; % set as 100 for percent


sample_tab = table(...
    [1:12]',...
    [datetime(1996,3,31), datetime(1996,3,31), datetime(2009,1,31),...
    datetime(1996,3,31), datetime(2000,1,31), datetime(2000,1,31),...
    datetime(1996,3,31), datetime(2000,1,31), datetime(1980,1,31),...
    datetime(1980,1,31), datetime(1988,2,28), datetime(1988,2,28),...
     ]',...
    [datetime(2019,12,31), datetime(2008,12,31), datetime(2019,12,31),...
    datetime(2024,4,30), datetime(2024,4,30), datetime(2019,12,31),...
    datetime(2024,1,31), datetime(2024,1,31), datetime(2024,1,31),...
    datetime(1996,3,31), datetime(1996,3,31), datetime(2024,5,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});


%{
sample_tab = table(...
    [1:3]',...
    [datetime(1986,1,31), datetime(1986,1,31), datetime(2009,1,31),...
     ]',...
    [datetime(2019,12,31), datetime(2000,12,31), datetime(2019,12,31),...
    ]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}


%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});
tab = table;
VariableNames = {'min_date', 'mid_date', 'max_date',...
    'alpha', 't_stat_alpha',...
    'beta', 't_stat_beta',...
    'corr', 'Chow'};


RowNames = arrayfun(@(lag) sprintf('lag%d', lag),...
    0:4, 'UniformOutput', false);


series_CPI = 'CP_NSA.M';

series_Rf = 'TSB_BAGR_MAKAM_01M.M';
series_tase_all_share = 'SPASTT01ILM661N'; % tase-all-share
series_tase_125 = 'MD_M137.M';
series_sm = series_tase_125;


%{
if ~vintage
    lib_ns_vintage = lib_ns;
end
%}

lib_EITB = fullfile(root, 'EITB_il');
lib_UITB = fullfile(root, 'UITB_il');
lib_tase_all_share = fullfile(lib_data, 'stock_markets');
lib_tase_125 = fullfile(lib_israel, 'SM', 'M', 'with_metadata');
lib_ir = fullfile(lib_israel, 'TSB_BAGR_MAKAM', 'M', 'with_metadata');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


data = struct(...
    'series', {...
    'CP_NSA.M',... CPI, not seasonally adjusted, BLS
    'CP_SA.M',... Price index of all urban consumers
    series_tase_all_share,...
    series_tase_125,...
    series_Rf,...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_cpi,...
    lib_tase_all_share,...
    lib_tase_125,...
    lib_ir,...
    %'',...
    },...
    'src', {...
    'FAME',...
    'FAME',...
    'FRED',...
    'FAME',...
    'FAME',...
    %'',...
    },...
    'VariableNames', cell(1,5)...
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
        %metadata = [metadata; metadata_];
        
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

% seasonally adjust CP_NSA for period not available from CBS, adding 12
% months for sa-estimation
CP_SA = sa_adj(T.CP_NSA0x2EM(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')+12), 12);

% set CP_SA as the official CBS series of the CPI
T.CP_SA = T.CP_SA0x2EM;

% complete this series with the above CP_SA
T.CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1) = ...
    CP_SA(1:find(~isnan(T.CP_SA0x2EM), 1, 'first')-1);

%{
% check sa_adj valid, pre 1996 data.
clf
hold on
plot(T.date, T.CP_SA)
plot(T.date, T.CP_SA0x2EM)
%plot(T.date, T.CP_NSA0x2EM)
%}



% Sample period
%for period = sample_tab.period(end-2)'
for period = sample_tab.period(end)'

    min_date = sample_tab.min_date(period);
    max_date = sample_tab.max_date(period);
        
    % ***************************************
    % TASE value weighted inedx

    idx = ismember({data.series}, series_sm);
    TASE = T{:, data(idx).VariableNames};
    
    % ***************************************
    % Risk free rate

    idx = ismember({data.series}, series_Rf);
    Rf = T{:, data(idx).VariableNames};

    % ***************************************
    % CPI

    CPI = T.CP_SA;
    T.Pi = [NaN;...
        CPI(2:end, :)./...
        CPI(1:end-1, :)];
    
    
    T.dTASE = [NaN;...
        TASE(2:end, :)./...
        TASE(1:end-1, :)];
    
   
    if log_approx
        
        %Pi = log(Pi)*pcent_mult;
        T.dTASE = log(T.dTASE)*pcent_mult;

    else

        %Pi = (Pi - 1)*pcent_mult;
        T.dTASE = (T.dTASE- 1)*pcent_mult;

    end
    
    T.Pi = (T.Pi - 1)*pcent_mult;
    T.RS = T.dTASE-T.Pi; % real stock return
    
    idx_sample = true(size(T, 1), 1);
    idx_sample = idx_sample & T.date >= min_date;
    idx_sample = idx_sample & T.date <= max_date;
    idx_sample = idx_sample & ~ismember(T.date, [...
        datetime(2008, 9, 30),...
        datetime(2008, 10, 31),...
        ]);
    
    T = T(idx_sample, :);

    mov_win = 60+1;
    for t = mov_win:size(T)

        min_date = T.date(t-mov_win+1);
        max_date = T.date(t);
        mid_date = T.date(t+1-(1+(mov_win -1)/2));
        
        idx_sample = T.date >= min_date & T.date <= max_date;
        idx_sample_pre = T.date >= min_date & T.date < mid_date;
        idx_sample_post = T.date > mid_date & T.date <= max_date;


        Y = T.dTASE(idx_sample);
        X = T.Pi(idx_sample);
        Y_pre = T.dTASE(idx_sample_pre);
        X_pre = T.Pi(idx_sample_pre);
        Y_post = T.dTASE(idx_sample_post);
        X_post = T.Pi(idx_sample_post);
    

        RowNames = {'Constant', 'Inflation'};
        X = [ones(size(X, 1), 1), X];

        N = size(X, 1);
        K = size(X, 2);
        
        N_pre = size(X_pre, 1)
        N_post = size(X_post, 1)
    
    
        % Full
        [B,BINT,R,~,STATS] = regress(Y, X);
        [B_pre, BINT_pre,R_pre,~, STATS_pre] = regress(Y_pre, X_pre);
        [B_post, BINT_post,R_post,~, STATS_post] = regress(Y_post, X_post);
    
        [B_, STDB_] = lscov(X, Y);
        [B_pre_, STDB_pre_] = lscov(X_pre, Y_pre);
        [B_post_, STDB_post_] = lscov(X_post, Y_post);

        SSR = sum(R.^2);
        SSR_pre = sum(R_pre.^2);
        SSR_post = sum(R_post.^2);

        C = (SSR - SSR_pre - SSR_post)/K/...
            ((SSR_pre + SSR_post)/(N_pre + N_post - 2*K));

        tab_row = table(...
            min_date, mid_date, max_date,...
            B(1), B_(1)./STDB_(1),...
            B(end), B_(end)./STDB_(end),...
            corr(X(:, end), Y), C,...
            'VariableNames', VariableNames);
        tab = [tab; tab_row];

        if mid_date == datetime(2008, 10, 31)
            [h, p] = jbtest(R)
            [h, p] = jbtest(R_pre)
            [h, p] = jbtest(R_post)
            clf
            hold on
            histogram(R_pre)
            histogram(R_post)
            
        end

    end

end

tab

%%
clf
set(gcf, 'Color', 'w')
set(gca, 'Fontsize', 14)
set(gcf, 'Position', [71   189   915   683])

%plot(tab.max_date, tab.beta)
%plot(tab.max_date, tab.alpha)
plot(tab.mid_date, tab.Chow)

set(gca, 'FontSize', 14)

%vline(datetime(2008, 10, 31), 'k--')
fname_fig = sprintf('chow_il_%dm_movwin.png', mov_win);
fpath_fig = fullfile(root, 'gfx', fname_fig);
export_fig(fpath_fig)






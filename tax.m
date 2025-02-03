%{

pp. 549
(6) Pi_t = -b_0 - b_1 dA_t - b_2 dRf_t + b_3 dM_t + eta_t

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
lib_ms = fullfile(lib_data, 'money_supply_usa');
lib_gnp = fullfile(lib_data, 'gnp');  % gross national product
lib_ea = fullfile(lib_data, 'economic_activity');  % economic activity
lib_cx = fullfile(lib_data, 'cx');  % capital expenditure
lib_ns = fullfile(lib_data, 'ns', 'level');  % capital expenditure
lib_tax = fullfile(lib_data, 'tax');  % capital expenditure


lib_EITB = fullfile(root, 'EITB');
lib_UITB = fullfile(root, 'UITB');

lib_ir = fullfile(lib_data, 'ir');
lib_pop = fullfile(lib_data, 'population');
lib_consumption = fullfile(lib_data, 'consumption_usa');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


% flag to override period in downstream scripts: false by default
lag_inflation = 0;
%log_approx = true;
log_approx = false;
freq = 'Q';

%{
% original dates
sample_tab = table(...
    [1:2]',...
    [datetime(1953,1,31), datetime(1954,1,31)]',...
    [datetime(1976,12,31), datetime(1976,2,28)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}
%
% extended dates
sample_tab = table(...
    [1:7]',...
    [datetime(1953,1,31), datetime(1953,1,31), datetime(2000,1,31),...
     datetime(1953,1,31), datetime(2007,7,31), datetime(1953,1,31),...
     datetime(1954,1,31)]',...
    [datetime(2019,12,31), datetime(1999,12,31), datetime(2019,12,31),...
    datetime(2007,7,31), datetime(2019,12,31), datetime(1977,12,31),...
    datetime(1976,12,31)]',...
    'VariableNames', {'period', 'min_date', 'max_date'});
%}

%sample_tab = sortrows(sample_tab, {'min_date', 'max_date'});

tab2 = [];


series_CPI = 'CPIAUCSL';
series_Rf = 'RF_FF';
%series_ms = 'BOGMBASE';
series_ms = {'CURRCIR', 'RESBALNS'};
series_ea = 'INDPRO';

% note: there are discrepancies between NIPA and FRED. use NIPA
series_cx_fred = 'BOGZ1FA105050005Q';
series_cx_nipa = 'BOGZ1FA105050005A_NIPA';

%series_gnp = 'GNPC96';
series_gnp = 'A791RX0Q048SBEA';

series_ns_fred = {'ESABSNNCB', 'RCSNNWMVBSNNCB', 'RCVSRNWMVBSNNCB'};
series_C = {'PCENDG', 'PCES'};
series_PD = {'PCEPINDG', 'PCEPIS'};


d = dir(lib_tax);
d = d(~[d.isdir]);
d = d(cellfun(@(c) contains(c, '.csv'), {d.name}));

% load data to T_ and join it 
for i = 1:numel(d)
    
    series = strrep(d(i).name, '.csv', '');
    
    [T_, metadata_] = load_data_files(...
        series, lib_tax, 'FRED');
    
    % robust date, start of month
    aux_date = T_.date+1;
    if all(aux_date.Day == 1)
        %T_.date = aux_date-calmonths(1);
        T_.date = aux_date;
    end
    if strcmp(series, series_CPI)
        % Boons lag inflation one period, but Fama does not.
        T_.date = T_.date + calmonths(lag_inflation);
    end
    
    if ~exist('T', 'var')
        T = T_;
        metadata = metadata_;
    else
        T = outerjoin(T, T_, 'Keys', {'date'},...
            'MergeKeys', true, 'Type', 'left');
        metadata = [metadata; metadata_];
    end
    

end

if strcmp(freq, 'Q')
    %pass
elseif strcmp(freq, 'M')
    T = table2timetable(T);
    T = retime(T, 'monthly', 'previous');
    T = timetable2table(T);
elseif strcmp(freq, 'A')
    T = table2timetable(T);
    T = retime(T, 'yearly', 'lastvalue');
    T = timetable2table(T);
else
    error('freq must be M, Q or A')
end

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month

%T = fillmissing(T, 'previous');


%%
clc

metadata = sortrows(metadata, 'series');
metadata

% TXPLBSNNCB, BOGZ1FL103178005Q, Level
% NCBPALQ027S Transactions
% do not mix level and transaction data. one is flow, the other is stock.
grp1 = {'BOGZ1FL103178005Q', 'TXPLBSNNCB'};
grp1a = {'BOGZ1FL103178005Q', 'BOGZ1FL143178005Q'};

% BOGZ1FL103178005Q, total taxes, liability
% BOGZ1FL143178005Q, taxes payable

%grp2: profits after tax, roughly the same from different sources
grp2 = {'NFCPATAX', 'A466RC1A027NBEA'};
grp3 = {'A464RC1Q027SBEA', 'A466RC1A027NBEA'};
grp4 = {'B3190C0A144NBEA', 'J3151C0A144NBEA', 'B3151C0A144NBEA'};
grp5 = {'A054RC1A027NBEA', 'BOGZ1FL143178005Q'};
grp6 = {'BOGZ1FL143178005Q', 'BOGZ1FL103178005Q'};


grp = sort(grp6);
% before 2015, they are the same.

%grp = sort(grp3);

clf
hold on
plot(T.date, T{:, grp})
%plot(T.date, sum(T{:, {'N3152C0A144NBEA', 'N3158C0A144NBEA'}}, 2))
%{
a=arrayfun(@(t) prod(T{2:t, grp}./T{1:t-1, grp}, 'omitnan')', 10:size(T,1),...
    'UniformOutput', false);
a = cell2mat(a);
b = T{:, grp(2)}./T{:, grp(1)};
%}
%plot(T.date(2:end), (T{2:end, grp}./T{1:end-1, grp}-1)*100)
%plot(T.date(2:end), T{2:end, grp}./T{1:end-1, grp})
%plot(T.date(2:end), prod(T{2:end, grp}./T{1:end-1, grp}))
%plot(T.date(2:end), cumsum((T{2:end, grp}./T{1:end-1, grp}-1)*100, 'omitnan'))
%plot(T.date(10:end), a')
%plot(T.date, b)

legend(metadata.desc{ismember(metadata.series, grp)})


%{
In 1998, NIPAs are based on NAICS, production oriented, compared to
previous classification using sic (demand oriented).

I've made a little exercise to see which series are direct continuation of
others, and which are not.

It is difficult to construct a consistent series for: ``Financial and
real-estate sector''. There seems to be a big gap, even when adding up the
real-estate industries in the period 1998-2021.

Resolution: use the series in grp4, and extend them using growth rates in
sum(T{:, {'N3152C0A144NBEA', 'N3158C0A144NBEA'}}, 2)

The same can be done for other series.

Compare this to the FRED non-fin series



%}

return

a = T{:, grp4};
a = sum(a, 2, 'omitnan')./sum(~isnan(a), 2);

b = sum(T{:, {'N3152C0A144NBEA', 'N3158C0A144NBEA'}}, 2, 'omitnan');
a(a==0) = NaN;
b(b==0) = NaN;

da = a(2:end)./a(1:end-1);
db = b(2:end)./b(1:end-1);

db(isnan(db)) = da(isnan(db));

fin = [a(1); arrayfun(@(i) a(1)*prod(db(1:i)), 1:size(db,1))'];
corp = T.A054RC1A027NBEA;
nonfin1 = corp-fin;
nonfin2 = T.BOGZ1FL143178005Q;% FRED series
nonfin3 = T.A464RC1Q027SBEA-T.A466RC1A027NBEA; % pbt-pat=tax


clf
hold on
plot(T.date, nonfin1)
plot(T.date, nonfin2)
plot(T.date, nonfin3)



return


%%

% compare grp3(1)-grp3(2) with grp1(2)
% they do not match, TXPLBSNNCB is much smaller. it's only other liab..

clf
hold on

plot(T.date, T{:, 'A464RC1Q027SBEA'}-T{:, 'A466RC1A027NBEA'})
plot(T.date, T{:, 'BOGZ1FL143178005Q'})
% these are close, but not enough. 


%%



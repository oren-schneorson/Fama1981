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
lib_profits = fullfile(lib_data, 'profits');  % capital expenditure


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


d = dir(lib_profits);
d = d(~[d.isdir]);
d = d(cellfun(@(c) contains(c, '.csv'), {d.name}));

% load data to T_ and join it 
for i = 1:numel(d)
    
    series = strrep(d(i).name, '.csv', '');
    
    [T_, metadata_] = load_data_files(...
        series, lib_profits, 'FRED');
    
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
metadata(:, {'series', 'desc', 'freq', 'sa'})

% BOGZ1FA146060035Q, includes something in addition?

grp1 = {'BOGZ1FA146060035Q', 'BOGZ1FU106060005Q', 'BOGZ1FU106060035Q'};
grp2 = {'BOGZ1FU106060005Q', 'BOGZ1FU106060035Q'};
grp3 = {'BOGZ1FU106060005Q', 'BOGZ1FU106060035Q', 'NCBIVDQ027S'};
grp4 = {'BOGZ1FU106060005Q', 'NCBIVDQ027S'};
grp5 = {'A464RC1Q027SBEA', 'BOGZ1FU106060005Q'};
grp6 = {'A464RC1Q027SBEA', 'NFCPATAX'};
grp7 = {'A053RC1A027NBEA', 'BOGZ1FA796060005A'};
grp8 = {'A3051C0A144NBEA', 'J3051C0A144NBEA', 'A3050C0A144NBEA'};
grp9 = {'N3093C0A144NBEA', 'B3053C0A144NBEA', 'J3053C0A144NBEA'}; % federal reserve banks
grp10 = {'N3061C0A144NBEA', 'B3059C0A144NBEA', 'J3059C0A144NBEA', 'A3059C0A144NBEA'}; % real estate
grp11 = {'N3057C0A144NBEA', 'B3056C0A144NBEA', 'J3056C0A144NBEA', 'A3053C0A144NBEA'}; % security and commodity brokers
grp12 = {'B3091C0A144NBEA', 'J3052C0A144NBEA', 'A3052C0A144NBEA'}; % banking / depository institutions
grp13 = {'N3058C0A144NBEA', 'B3057C0A144NBEA', 'B3058C0A144NBEA'}; % 
grp14 = {'N3062C0A144NBEA', 'B3060C0A144NBEA'}; % 
grp15 = {'N3062C0A144NBEA', 'A464RC1Q027SBEA'}; % 


% N3060C0A144NBEA + N3092C0A144NBEA

grp = sort(grp15);
% before 2015, they are the same.

clf
hold on
%plot(T.date, [T{:, 'N3092C0A144NBEA'}, [T{:, 'A3051C0A144NBEA'}-sum(T{:, {'B3059C0A144NBEA', 'B3060C0A144NBEA'}}, 2)]])
%plot(T.date, [sum(T{:, {'N3060C0A144NBEA', 'N3092C0A144NBEA'}}, 2), [T{:, 'A3051C0A144NBEA'}-T{:, 'B3059C0A144NBEA'}]])
%plot(T.date, [sum(T{:, {'N3061C0A144NBEA', 'N3092C0A144NBEA'}}, 2), T{:, grp8}])
%plot(T.date, [sum(T{:, {'N3093C0A144NBEA', 'N3056C0A144NBEA'}}, 2), T{:, grp12}]) % not quite
%{
% insurance.
plot(T.date, [T{:, 'N3058C0A144NBEA'},...
    sum(T{:, {'B3057C0A144NBEA', 'B3058C0A144NBEA'}}, 2),...
    sum(T{:, {'J3057C0A144NBEA', 'J3058C0A144NBEA'}}, 2),...
    sum(T{:, {'A3055C0A144NBEA', 'A3056C0A144NBEA'}}, 2),...
    ])
%}
plot(T.date, T{:, grp})


legend(metadata.desc{ismember(metadata.series, grp)})

%{
In 1998, NIPAs are based on NAICS, production oriented, compared to
previous classification using sic (demand oriented).

I've made a little exercise to see which series are direct continuation of
others, and which are not.

It is difficult to construct a consistent series for: ``Financial and
real-estate sector''. There seems to be a big gap, even when adding up the
real-estate industries in the period 1998-2021.

Resolution: use the series in grp8, and extend them using growth rates in
sum(T{:, {'N3061C0A144NBEA', 'N3092C0A144NBEA'}}, 2)

The same can be done for other series



%}


a = T{:, grp8};
a = sum(a, 2, 'omitnan')./sum(~isnan(a), 2);

b = sum(T{:, {'N3061C0A144NBEA', 'N3092C0A144NBEA'}}, 2, 'omitnan');
a(a==0) = NaN;
b(b==0) = NaN;

da = a(2:end)./a(1:end-1);
db = b(2:end)./b(1:end-1);

db(isnan(db)) = da(isnan(db));

fin = [a(1); arrayfun(@(i) a(1)*prod(db(1:i)), 1:size(db,1))'];
corp = T.A053RC1A027NBEA;
nonfin1 = corp-fin;

nonfin2 = T.BOGZ1FA106060005Q;
nonfin3 = T.BOGZ1FA106060035Q;
nonfin4 = T.A464RC1A027NBEA;
nonfin5 = T.A464RC1Q027SBEA;

clf
hold on
plot(nonfin1)
plot(nonfin5)
return

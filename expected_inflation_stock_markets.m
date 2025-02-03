
%{
Estimating how expected inflation and excess returns co-move.
U.S. data.

Stock market is taken to be S&P-500 from Yahoo Finance figures.
Risk free rate is 3 months U.S. T-bill.
Expected inflation is taken from TIPS, which are in daily
frequency, so that all other data is kept to daily.

Last check: 2023-01-17



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


series_SM = 'SP500';
fname = [series_SM, '.xlsx'];
fpath = fullfile(lib_data, fname);

T_ = readtable(fpath);
T_.Properties.VariableNames{1} = 'date';
T_.Properties.VariableNames{...
    ismember(T_.Properties.VariableNames, 'AdjClose')} = series_SM;
%T_.ret = (T_.Close./T_.Open-1)*100;
T_ = sortrows(T_, 'date');
T_.ret = [NaN; T_.(series_SM)(2:end)./T_.(series_SM)(1:end-1)];
T_.ret = (T_.ret-1)*100;

T = T_;
clear T_

series_RF = 'DTB3';
fname = [series_RF, '.xlsx'];
fpath = fullfile(lib_data, 'ir', 'DTB', fname);
T_ = readtable(fpath);
T_.Properties.VariableNames = {'date', series_RF};

%T = outerjoin(T, T_, 'Keys', {'date'}, 'MergeKeys', true);
T = innerjoin(T, T_, 'Keys', {'date'});
clear T_

series = 'CPI';
fname = [series, '.xlsx'];
fpath = fullfile(lib_data, fname);

T_ = readtable(fpath);
T_.observation_date = datetime(T_.observation_date);
T_.Properties.VariableNames{1} = 'date';
series_CPI = 'USACPIALLMINMEI';
keep = {'date', series_CPI};
T_ = T_(:, keep);
T_.date = T_.date-1;


T = outerjoin(T, T_, 'Keys', {'date'}, 'MergeKeys', true);
clear T_

x = T{:, series_CPI};
nanx = isnan(x);
t = 1:numel(x);
x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
T{:, series_CPI} = x;

%% compute 1M, 1Y and 5Y actual inflation
clc
%{
aux = T.date.Month-1;


[datetime(...
    T.date.Year -(aux==0),...
    T.date.Month-1+(aux==0)*12,...
    T.date.Day)]
%}

% for each date in T, find the closest date 1 month before
aux = closest_before_date(T.date, 1, 'month');
aux = T.date' == aux;
NaNs = sum(aux, 2) == 0;
[~, aux] = find(aux);
aux_ = NaN(size(T,1), 1);
aux_(~NaNs) = T{~NaNs, series_CPI}./T{aux, series_CPI};
T.Pi_1M = aux_;

% for each date in T, find the closest date 1 year before
aux = closest_before_date(T.date, 1, 'year');
aux = T.date' == aux;
NaNs = sum(aux, 2) == 0;
[~, aux] = find(aux);
aux_ = NaN(size(T,1), 1);
aux_(~NaNs) = T{~NaNs, series_CPI}./T{aux, series_CPI};
T.Pi_1Y = aux_;

% for each date in T, find the closest date 5 years before
aux = closest_before_date(T.date, 5, 'year');
aux = T.date' == aux;
NaNs = sum(aux, 2) == 0;
[~, aux] = find(aux);
aux_ = NaN(size(T,1), 1);
aux_(~NaNs) = T{~NaNs, series_CPI}./T{aux, series_CPI};
T.Pi_5Y = aux_;
T.Pi_5Y = T.Pi_5Y.^(1/5);



T.Pi_1M = (T.Pi_1M-1)*100;
T.Pi_1Y = (T.Pi_1Y-1)*100;
T.Pi_5Y = (T.Pi_5Y-1)*100;

T = T(~isnan(T.(series_CPI)),:); % gets rid of obs in 2022, no data on CPI
T.([series_SM, '_real']) = T.(series_SM)./T.(series_CPI);
T.ret_real = [NaN;...
    T.([series_SM, '_real'])(2:end)./...
    T.([series_SM, '_real'])(1:end-1)];
T.ret_real = (T.ret_real-1)*100;

fname = 'expected_inflation_TIPS.xlsx';
fpath = fullfile(lib_data, fname);


T_ = readtable(fpath);
T_.Properties.VariableNames{1} = 'date';
InputFormat = 'yyyy-MM-dd';
T_.date = datetime(T_.date, 'InputFormat', InputFormat);
T_ = sortrows(T_, 'date');

for v = T_.Properties.VariableNames(2:end)
   
    v = char(v);
    v_new = ['D', v];
    T_(:, v_new) = table([NaN; T_{2:end, v}-T_{1:end-1, v}]);
    
end

T = innerjoin(T, T_, 'Keys', {'date'});
clear T_

%%
T.DTB3_real = T.DTB3 - ((T.T5YIE/100+1).^(1/4)-1)*100;

% from annual to daily rate
T.DTB3 = 1 + T.DTB3/100;
T.DTB3 = T.DTB3.^(1/252);
T.DTB3 = (T.DTB3-1)*100;

T.DTB3_real = 1 + T.DTB3_real/100;
T.DTB3_real = T.DTB3_real.^(1/252);
T.DTB3_real = (T.DTB3_real-1)*100;


T.xret = T.ret-T.DTB3;
%T.xret_real = T.ret_real-T.DTB3; % minus nominal rate
T.xret_real = T.ret_real-T.DTB3_real; %minus real rate


return

%%
clc
Y_var = {'ret_real'};
%Y_var = {'AdjClose'};
X_var = {'DT5YIE','DT5YIFR'};
%X_var = {'T5YIE','T5YIFR','DT5YIE','DT5YIFR'};
%X_var = {'T5YIE','T5YIFR'};
%X_var = {'DT5YIE','DT5YIFR','DT10YIE'};
%X_var = {'T5YIE','T5YIFR','T10YIE'};
%X_var = {'T5YIE','T5YIFR','T10YIE','DT5YIE','DT5YIFR','DT10YIE'};

date_start = datetime(2011,1,1);
date_end = datetime(2019,12,31);
%date_start = datetime(2000,1,1);
%date_end = datetime(2006,12,31);
%date_start = min(T.date);
%date_end = max(T.date);

idx = T.date >= date_start & T.date <= date_end;

Y = T{idx, Y_var};
X = T{idx, X_var};
X = [ones(size(X,1), 1), X];

idx = all(~isnan([Y,X]), 2);
X = X(idx,:);
Y = Y(idx,:);

[B,BINT,R,RINT,STATS] = regress(Y,X);
[BINT(:,1), B, BINT(:,2)]

%% moving window regression
clc

constant = true;
%constant = false;
m_window = 60;
y_window = m_window/12;
d_window = round(365*y_window);

Y_var = {'ret_real'};
%Y_var = {'AdjClose'};
%X_var = {'DT5YIE','DT5YIFR'};
X_var = {'DT5YIE'};
%X_var = {'T5YIE','T5YIFR','DT5YIE','DT5YIFR'};
%X_var = {'T5YIE','T5YIFR'};
%X_var = {'DT5YIE','DT5YIFR','DT10YIE'};
%X_var = {'T5YIE','T5YIFR','T10YIE'};
%X_var = {'T5YIE','T5YIFR','T10YIE','DT5YIE','DT5YIFR','DT10YIE'};

results = struct;
results.DATE = NaT;
results.B = NaN(2,1);
results.BINT1 = NaN(2,1);
results.BINT2 = NaN(2,1);
results.R = NaN;
%results.RINT = NaN;
results.STATS = NaN(4,1);

counter = 0;
for d = T.date(1)+d_window:1:T.date(end)
    counter = counter+1;
    date_start = d-d_window;
    date_end = d;

    idx = T.date >= date_start & T.date <= date_end;

    Y = T{idx, Y_var};
    X = T{idx, X_var};
    if constant
    	X = [ones(size(X,1), 1), X];
    end

    idx = all(~isnan([Y,X]), 2);
    X = X(idx,:);
    Y = Y(idx,:);

    [B,BINT,R,~,STATS] = regress(Y,X);
    results(counter).DATE = d;
    results(counter).B = B;
    results(counter).BINT1 = BINT(:,1);
    results(counter).BINT2 = BINT(:,2);
    results(counter).R = R;
    %results(counter).RINT = RINT;
    results(counter).STATS = STATS';
    %[BINT(:,1), B, BINT(:,2)]
    

end

%
% plot results
clc
clf
FaceAlpha = .4;
lehman = datetime(2008,9,15);

coef = 1+constant;

hold on
DATE = [results.DATE];
B = [results.B];
pl = plot(DATE, B(coef,:));

BINT1 = [results.BINT1];
BINT2 = [results.BINT2];
BINT = [BINT2(coef,:); BINT1(coef,:)-BINT2(coef,:)];

%BINT = flip(BINT);
BINT = BINT';
%plot(DATE, BINT1(coef,:))
ar = area(BINT, 'LineStyle', 'none');
ar(1).FaceColor = 'none';
%ar(2).FaceColor = 'b';
ar(2).FaceColor = pl.Color;
ar(2).FaceAlpha = FaceAlpha;
vline(lehman)

%ar.LineStyle = 'none';
%colors

leg = [X_var, {'', '95%-Confidence interval'}];
%leg = [{'Constant'}, X_var];

legend(leg, 'Location', 'North')
tit = [char(Y_var), '_', char(X_var)];
tit = [char(tit), sprintf('; %d-month moving window daily', m_window)];

fname = [tit, '.png'];
fpath = fullfile(root, 'gfx', fname);
tit = strrep(tit, '_', '\_');
%title(tit)

export_fig(fpath)
return
%% plot series

Y_var = {'AdjClose_real'};
X_var = {'T5YIE','T5YIFR','T10YIE'};


date_start = min(T.date);
date_end = max(T.date);

idx = T.date >= date_start & T.date <= date_end;
Y = T{idx, Y_var};
X = T{idx, X_var};
date = T.date(idx);

clf
hold on
plot(date, Y)
yyaxis right
plot(date, X)

legend([Y_var, X_var])

%% plot inflation and expected inflation

clc
clf

var_Pi = 'Pi_5Y';
var_ex_Pi = 'T5YIE';
YLim = [-3 10];

idx = ~isnat(T.date);

hold on
plot(T.date(idx), T.(var_Pi)(idx))
ylim(YLim)
ylabel('Percent')
legend('Location', 'Northwest')

plot(T.date(idx), T.(var_ex_Pi)(idx))
ylabel('Percent')
ylim(YLim)
legend({...
    'Inflation (CPI, 1 year backward)',...
    'Expected inflation (5 year forward)'})

%title({'Inflation and expected inflation'})
ax = gca;
ax.FontSize = 14;
%ax.XAxis.TickLabelRotation = 30;

fname = 'Inflation (1 year backward) and Expected inflation (5 years forward)_USA.png';
fpath = fullfile(root, 'gfx', fname);

xlim([datetime(1998,9,2), datetime(2022,1,31)])
ax.XTick = [...
    datetime(1998,9,2),...
    datetime(2004,2,12),...
    datetime(2009,8,15),...
    datetime(2015,2,5),...
    datetime(2020,7,28)];

ax.XTickLabel = {'1998','2004','2009','2015','2020'}

set(gcf,'color','w');
export_fig(fpath)


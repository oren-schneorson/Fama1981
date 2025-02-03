%{
Estimating how expected inflation and excess returns co-move.
U.S. data.

Stock market is taken from Shiller data.
Risk free rate is 3 months U.S. T-bill.
Expected inflation is taken from the Cleveland model, which are in monthly
frequency, so that all other data is transformed to monthly.

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

% ie_data_with_TRCAPEseries = 'ie_data_with_TRCAPE';
fname = 'ie_data_with_TRCAPE.xlsx';
fpath = fullfile(lib_data, 'stock_markets', fname);

Range = 'A2:N1775';
Sheet = 'Data';

T_ = readtable(fpath, 'Sheet', Sheet, 'Range', Range);
T_.Properties.VariableNames{1} = 'date';
T_.date = num2str(T_.date,'%.2f');
T_.date = cellstr(num2str(T_.date));


InputFormat = 'yyyy.MM';
T_.date = datetime(T_.date, 'InputFormat', InputFormat);
T_ = sortrows(T_, 'date');

series_SM = 'P';
T_.ret = [NaN; T_.(series_SM)(2:end)./T_.(series_SM)(1:end-1)];
T_.ret = (T_.ret-1)*100;

T = T_;
clear T_

series = 'DTB3';
fname = [series, '.xlsx'];
fpath = fullfile(lib_data, 'ir', 'DTB', fname);
T_ = readtable(fpath);
T_.Properties.VariableNames = {'date', series};

T_ = table2timetable(T_);
T_ = retime(T_, 'daily', 'previous');
T_ = retime(T_, 'monthly');
T_ = timetable2table(T_);

T = outerjoin(T, T_, 'Keys', {'date'}, 'MergeKeys', true);
clear T_



var_CPI = 'CPI';
keep = {'date', var_CPI};

T.ret_real = [NaN; T.RTRP(2:end)./T.RTRP(1:end-1)];
T.ret_real = (T.ret_real-1)*100;


fname = 'expected_inflation_Cleveland_model.xlsx';
fpath = fullfile(lib_data, fname);


T_ = readtable(fpath);
T_.observation_date = datetime(T_.observation_date);
T_.Properties.VariableNames{1} = 'date';
T_ = sortrows(T_, 'date');

for v = T_.Properties.VariableNames(2:end)
   
    v = char(v);
    v_new = ['D', v];
    T_(:, v_new) = table([NaN; T_{2:end, v}-T_{1:end-1, v}]);
    
end

T = innerjoin(T, T_, 'Keys', {'date'});
clear T_


% subtract expected inflation for the next quarter
T.DTB3_real = T.DTB3 - ((T.EXPINF1YR/100+1).^(1/4)-1)*100;

% from annual to monthly rate
T.DTB3 = 1 + T.DTB3/100;
T.DTB3 = T.DTB3.^(1/12);
T.DTB3 = (T.DTB3-1)*100;

T.DTB3_real = 1 + T.DTB3_real/100;
T.DTB3_real = T.DTB3_real.^(1/12);
T.DTB3_real = (T.DTB3_real-1)*100;


T.xret = T.ret-T.DTB3;
T.xret_real = T.ret_real-T.DTB3; % minus nominal rate
%T.xret_real = T.ret_real-T.DTB3_real; minus real rate


%%
clc
constant = 1;
Y_var = {'ret'};
%Y_var = {'AdjClose_real'};
X_var = {'DEXPINF3YR'};
%X_var = {'T5YIE','T5YIFR','DT5YIE','DT5YIFR'};
%X_var = {'T5YIE','T5YIFR'};
%X_var = {'DT5YIE','DT5YIFR','DT10YIE'};
%X_var = {'T5YIE','T5YIFR','T10YIE'};
%X_var = {'T5YIE','T5YIFR','T10YIE','DT5YIE','DT5YIFR','DT10YIE'};

%date_start = datetime(2011,1,1);
%date_end = datetime(2019,12,31);
%date_start = datetime(2000,1,1);
date_end = datetime(1995,12,31);
date_start = min(T.date);
%date_end = max(T.date);

idx = T.date >= date_start & T.date <= date_end;

Y = T{idx, Y_var};
X = T{idx, X_var};
if constant
    X = [ones(size(X,1), 1), X];
end

idx = all(~isnan([Y,X]), 2);
X = X(idx,:);
Y = Y(idx,:);

[B,BINT,R,RINT,STATS] = regress(Y,X);
[BINT(:,1), B, BINT(:,2)]
SS_res = sum(R.^2); % sum of squared residuals
SS_tot = sum((Y-mean(Y)).^2); % total sum of squares
R_sq = 1-SS_res/SS_tot
STATS(1)




%% moving window regression
clc

constant = true;
%constant = false;
y_window = 5;
m_window = y_window*12;

Y_var = {'xret'};
%Y_var = {'AdjClose'};
%X_var = {'DT5YIE','DT5YIFR'};
%X_var = {'DEXPINF5YR'};
X_var = {'DEXPINF10YR'};

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
for d = T.date(m_window:end)'
    if d >= datetime(2020,9,1)
        continue
    end
    counter = counter+1;
    date_start = datetime(d.Year-y_window, d.Month, d.Day);
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
ar = area(DATE, BINT, 'LineStyle', 'none');
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
tit = Y_var;
tit = [char(tit), sprintf('; %d-month moving window monthly', m_window)];

fname = [tit, '.png'];
fpath = fullfile(root, 'gfx', fname);
tit = strrep(tit, '_', '\_');

title(tit)

%export_fig(fpath)
return
%% plot series

clc
clf
Y_var = series_SM;
X_var = {'EXPINF5YR'};


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


leg = strrep([Y_var, X_var], '_', '\_');
legend(leg)










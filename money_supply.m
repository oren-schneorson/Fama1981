%{

Estimating the nominal-real covariance (NRC, Boons et al.)


Data:
1. Monthly inflation, all urban consumers
2. Monthly consumption, BEA
3. Price deflator, BEA
4. Population numbers, BEA

Boons et al. use ARMA(1,1) model for inflation to gauge inflation
innovations.

Innovations:
Look at the change in the 12-months forecast; 
(not the forecast error)
Also can look at the 1M forecast error.


References:

Boons, M., Duarte, F., De Roon, F., & Szymanowska, M. (2020).
Time-varying inflation risk and stock returns. 
Journal of Financial Economics, 136(2), 444-470.

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

lib_ir = fullfile(lib_data, 'ir');
lib_pop = fullfile(lib_data, 'population');
lib_consumption = fullfile(lib_data, 'consumption_usa');
lib_fama_french_factors = fullfile(lib_data, 'fama_french_factors');


series_CPI = 'CPIAUCSL';
series_Rf = 'RF_FF';
series_C = {'PCENDG', 'PCES'};
series_PD = {'PCEPINDG', 'PCEPIS'};


% src: BI, Bank of Israel old website; FAME
data = struct(...
    'series', {...
     series_CPI,... Price index of all urban consumers
    'BOGMBASE',...
    'BOGMBBM',...
    'BORROW',...
    'CURRCIR',...
    'M1SL',...
    'M2SL',...
    'MBCURRCIR',...
    'NONBORRES',...
    'TOTRESNS',...
    %'',...
    },...
    'lib_data', {...
    lib_cpi,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    lib_ms,...
    %'',...
    },...
    'src', {...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    'FRED',...
    %'',...
    },...
    'VariableNames', cell(1,10)...
    );


% declare containers
metadata = struct;
for src = unique({data.src})
    src = char(src);
    metadata.(src) = table;
end


% load data to T_ and join it 
for i = 1:size(data, 2)
    
    lib_data = data(i).lib_data;
    [T_, metadata_] = load_data_files(data(i).series, lib_data, data(i).src);
    data(i).VariableNames = T_.Properties.VariableNames(2);
    
    % robust date, start of month
    aux_date = T_.date+1;
    if all(aux_date.Day == 1)
        %T_.date = aux_date-calmonths(1);
        T_.date = aux_date;
    end
    if strcmp(data(i).series, series_CPI)
        % Boons lag inflation one period
        %T_.date = T_.date + calmonths(1);
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

%T.date = T.date + calmonths(1); % to end of month
T.date = T.date - caldays(1); % to end of month


return

%%

%{
BOGMBASE ~= BOGMBBM + CURRCIR
but they are close, up to 1% difference

%}
clf
plot(T.date, ((T.CURRCIR+T.BOGMBBM)./T.BOGMBASE-1)*100)



%%
%{
CURRCIR ~= MBCURRCIR
but they are close, up to 1% differnece
%}
clf
plot(T.date, (T.CURRCIR./T.MBCURRCIR-1)*100)




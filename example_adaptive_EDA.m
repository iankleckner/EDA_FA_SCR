% Ian Kleckner
% 2024/01/12

%% Setup
clc
clear all
clf

%% Import test dataset

data = readtable('data/example_data_Figure1.csv');
data_time_sec_Fig1 = data.Var1;
data_EDA_uS_Fig1 = data.Var2;

%% Test dataset
clf;
fprintf('\n-----------------------------------\nSTART\n');
RAP_test = 1.9;

SCR_Results = adaptiveEDA( data_time_sec_Fig1, data_EDA_uS_Fig1, RAP_test, true );

print(gcf, sprintf('out-Fig1_%0.2f.png', RAP_test), '-dpng')
fprintf('\nDone!\n')
fprintf('\nFound %d SCRs\n', length(SCR_Results.SCR_onset_time))

%% file example_L21.m
% this file shows the usage of Least_L21.m function 
% and study the group sparsity patterns. 
%
%% OBJECTIVE
% argmin_W { sum_i^t (0.5 * norm (Y{i} - X{i}' * W(:, i))^2) 
%            + rho1 * \|W\|_{2,1} + opts.rho_L2 * \|W\|_F^2}
%
%% LICENSE
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2011 - 2012 Jiayu Zhou and Jieping Ye 
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Jiayu Zhou via jiayu.zhou@asu.edu
%
%   Last modified on June 3, 2012.
%
%% Related papers
%
% [1] Evgeniou, A. and Pontil, M. Multi-task feature learning, NIPS 2007.
% [2] Liu, J. and Ye, J. Efficient L1/Lq Norm Regularization, Technical
% Report, 2010.
%
clear;
clc;
close;
addpath('C:/Users/hnr59/Documents/MATLAB/jiayuzhou-MALSAR-fb97515/MALSAR/functions/joint_feature_learning'); % load function
addpath('C:/Users/hnr59/Documents/MATLAB/jiayuzhou-MALSAR-fb97515/MALSAR/utils/'); % load utilities
%load('data/school.mat'); % load sample data.

%This data inputs the data subsets I created

M_L1 = readmatrix("df_computed_mental_health_L1_0days.csv");
M_L2 = readmatrix("df_computed_mental_health_L2_1to13days.csv");
M_L3 = readmatrix("df_computed_mental_health_L3_14plusdays.csv");

%S = cell2struct(num2cell(M, 1), compose('Var%d', 1:size(M, 2)), 2)
%Create cells for X and Y
%For this example, simply copy over the full table and multiply each value
%by 2 for X{2}. Also, for a test run, only include the first 149 rows
X = {M_L1(4:end, [4:end-1]), M_L2(4:end, [4:end-1]), M_L3(4:end, [4:end-1])}
Y = {M_L1(4:end, end), M_L2(4:end, end), M_L3(4:end, end)}
%X = M(:, [2:50 52:end]);
%Y = M(:, 51);
d = size(X{1}, 2);  % dimensionality.
lambda = [200 :300: 1500];
%rng('default');     % reset random generator. Available from Matlab 2011.
opts.init = 0;      % guess start point from data. 
opts.tFlag = 1;     % terminate after relative objective value does not changes much.
opts.tol = 10^-5;   % tolerance. 
opts.maxIter = 1000; % maximum iteration number of optimization.
sparsity = zeros(length(lambda), 1);
log_lam  = log(lambda);
for i = 1: length(lambda)
    [W funcVal] = Least_L21(X, Y, lambda(i), opts);
    % set the solution as the next initial point. 
    % this gives better efficiency. 
    opts.init = 1;
    opts.W0 = W;
    sparsity(i) = nnz(sum(W,2 )==0)/d;
end
% draw figure
h = figure;
plot(log_lam, sparsity);
xlabel('log(\rho_1)')
ylabel('Row Sparsity of Model (Percentage of All-Zero Columns)')
title('Row Sparsity of Predictive Model when Changing Regularization Parameter');
set(gca,'FontSize',12);
print('-dpdf', '-r100', 'LeastL21Exp');
%% 
W_toCsv = [(0:size(W,1)-1)' W];
writematrix(W_toCsv,'W_MTL_Least_L21__MENT14D.xlsx');

%%

% The data is set for _MENT14D as an output (last column) so the states
% data would need to be fix.
clear;
clc;
close;
addpath('C:/Users/hnr59/Documents/MATLAB/jiayuzhou-MALSAR-fb97515/MALSAR/functions/joint_feature_learning'); % load function
addpath('C:/Users/hnr59/Documents/MATLAB/jiayuzhou-MALSAR-fb97515/MALSAR/utils/'); % load utilities
addpath("C:/Users/hnr59/Documents/MATLAB/states_data/")
%load('data/school.mat'); % load sample data.


M_LS1 = readmatrix('df_state_L1.csv');
M_LS2 = readmatrix('df_state_L2.csv');
M_LS4 = readmatrix('df_state_L4.csv');
M_LS5 = readmatrix('df_state_L5.csv');
M_LS6 = readmatrix('df_state_L6.csv');
M_LS8 = readmatrix('df_state_L8.csv');
M_LS9 = readmatrix('df_state_L9.csv');
M_LS10 = readmatrix('df_state_L10.csv');
M_LS11 = readmatrix('df_state_L11.csv');
M_LS13 = readmatrix('df_state_L13.csv');
M_LS15 = readmatrix('df_state_L15.csv');
M_LS16 = readmatrix('df_state_L16.csv');
M_LS17 = readmatrix('df_state_L17.csv');
M_LS18 = readmatrix('df_state_L18.csv');
M_LS19 = readmatrix('df_state_L19.csv');
M_LS20 = readmatrix('df_state_L20.csv');
M_LS21 = readmatrix('df_state_L21.csv');
M_LS22 = readmatrix('df_state_L22.csv');
M_LS23 = readmatrix('df_state_L23.csv');
M_LS24 = readmatrix('df_state_L24.csv');
M_LS25 = readmatrix('df_state_L25.csv');
M_LS26 = readmatrix('df_state_L26.csv');
M_LS27 = readmatrix('df_state_L27.csv');
M_LS28 = readmatrix('df_state_L28.csv');
M_LS29 = readmatrix('df_state_L29.csv');
M_LS30 = readmatrix('df_state_L30.csv');
M_LS31 = readmatrix('df_state_L31.csv');
M_LS32 = readmatrix('df_state_L32.csv');
M_LS33 = readmatrix('df_state_L33.csv');
M_LS34 = readmatrix('df_state_L34.csv');
M_LS35 = readmatrix('df_state_L35.csv');
M_LS36 = readmatrix('df_state_L36.csv');
M_LS37 = readmatrix('df_state_L37.csv');
M_LS38 = readmatrix('df_state_L38.csv');
M_LS39 = readmatrix('df_state_L39.csv');
M_LS40 = readmatrix('df_state_L40.csv');
M_LS41 = readmatrix('df_state_L41.csv');
M_LS42 = readmatrix('df_state_L42.csv');
M_LS44 = readmatrix('df_state_L44.csv');
M_LS45 = readmatrix('df_state_L45.csv');
M_LS46 = readmatrix('df_state_L46.csv');
M_LS47 = readmatrix('df_state_L47.csv');
M_LS48 = readmatrix('df_state_L48.csv');
M_LS49 = readmatrix('df_state_L49.csv');
M_LS50 = readmatrix('df_state_L50.csv');
M_LS51 = readmatrix('df_state_L51.csv');
M_LS53 = readmatrix('df_state_L53.csv');
M_LS54 = readmatrix('df_state_L54.csv');
M_LS55 = readmatrix('df_state_L55.csv');
M_LS56 = readmatrix('df_state_L56.csv');
M_LS66 = readmatrix('df_state_L66.csv');
M_LS72 = readmatrix('df_state_L72.csv');
M_LS78 = readmatrix('df_state_L78.csv');

%S = cell2struct(num2cell(M, 1), compose('Var%d', 1:size(M, 2)), 2)
%Create cells for X and Y
%For this example, simply copy over the full table and multiply each value
%by 2 for X{2}. Also, for a test run, only include the first 149 rows
X = { M_LS1(2:end, [2:end-1]),
 M_LS2(2:end, [2:end-1]),
 M_LS4(2:end, [2:end-1]),
 M_LS5(2:end, [2:end-1]),
 M_LS6(2:end, [2:end-1]),
 M_LS8(2:end, [2:end-1]),
 M_LS9(2:end, [2:end-1]),
 M_LS10(2:end, [2:end-1]),
 M_LS11(2:end, [2:end-1]),
 M_LS13(2:end, [2:end-1]),
 M_LS15(2:end, [2:end-1]),
 M_LS16(2:end, [2:end-1]),
 M_LS17(2:end, [2:end-1]),
 M_LS18(2:end, [2:end-1]),
 M_LS19(2:end, [2:end-1]),
 M_LS20(2:end, [2:end-1]),
 M_LS21(2:end, [2:end-1]),
 M_LS22(2:end, [2:end-1]),
 M_LS23(2:end, [2:end-1]),
 M_LS24(2:end, [2:end-1]),
 M_LS25(2:end, [2:end-1]),
 M_LS26(2:end, [2:end-1]),
 M_LS27(2:end, [2:end-1]),
 M_LS28(2:end, [2:end-1]),
 M_LS29(2:end, [2:end-1]),
 M_LS30(2:end, [2:end-1]),
 M_LS31(2:end, [2:end-1]),
 M_LS32(2:end, [2:end-1]),
 M_LS33(2:end, [2:end-1]),
 M_LS34(2:end, [2:end-1]),
 M_LS35(2:end, [2:end-1]),
 M_LS36(2:end, [2:end-1]),
 M_LS37(2:end, [2:end-1]),
 M_LS38(2:end, [2:end-1]),
 M_LS39(2:end, [2:end-1]),
 M_LS40(2:end, [2:end-1]),
 M_LS41(2:end, [2:end-1]),
 M_LS42(2:end, [2:end-1]),
 M_LS44(2:end, [2:end-1]),
 M_LS45(2:end, [2:end-1]),
 M_LS46(2:end, [2:end-1]),
 M_LS47(2:end, [2:end-1]),
 M_LS48(2:end, [2:end-1]),
 M_LS49(2:end, [2:end-1]),
 M_LS50(2:end, [2:end-1]),
 M_LS51(2:end, [2:end-1]),
 M_LS53(2:end, [2:end-1]),
 M_LS54(2:end, [2:end-1]),
 M_LS55(2:end, [2:end-1]),
 M_LS56(2:end, [2:end-1]),
 M_LS66(2:end, [2:end-1]),
 M_LS72(2:end, [2:end-1]),
 M_LS78(2:end, [2:end-1])};

Y = { M_LS1(2:end, end),
 M_LS2(2:end, end),
 M_LS4(2:end, end),
 M_LS5(2:end, end),
 M_LS6(2:end, end),
 M_LS8(2:end, end),
 M_LS9(2:end, end),
 M_LS10(2:end, end),
 M_LS11(2:end, end),
 M_LS13(2:end, end),
 M_LS15(2:end, end),
 M_LS16(2:end, end),
 M_LS17(2:end, end),
 M_LS18(2:end, end),
 M_LS19(2:end, end),
 M_LS20(2:end, end),
 M_LS21(2:end, end),
 M_LS22(2:end, end),
 M_LS23(2:end, end),
 M_LS24(2:end, end),
 M_LS25(2:end, end),
 M_LS26(2:end, end),
 M_LS27(2:end, end),
 M_LS28(2:end, end),
 M_LS29(2:end, end),
 M_LS30(2:end, end),
 M_LS31(2:end, end),
 M_LS32(2:end, end),
 M_LS33(2:end, end),
 M_LS34(2:end, end),
 M_LS35(2:end, end),
 M_LS36(2:end, end),
 M_LS37(2:end, end),
 M_LS38(2:end, end),
 M_LS39(2:end, end),
 M_LS40(2:end, end),
 M_LS41(2:end, end),
 M_LS42(2:end, end),
 M_LS44(2:end, end),
 M_LS45(2:end, end),
 M_LS46(2:end, end),
 M_LS47(2:end, end),
 M_LS48(2:end, end),
 M_LS49(2:end, end),
 M_LS50(2:end, end),
 M_LS51(2:end, end),
 M_LS53(2:end, end),
 M_LS54(2:end, end),
 M_LS55(2:end, end),
 M_LS56(2:end, end),
 M_LS66(2:end, end),
 M_LS72(2:end, end),
 M_LS78(2:end, end)}
%X = M(:, [2:50 52:end]);
%Y = M(:, 51);
d = size(X{1}, 2);  % dimensionality.
lambda = [200 :300: 1500];
%rng('default');     % reset random generator. Available from Matlab 2011.
opts.init = 0;      % guess start point from data. 
opts.tFlag = 1;     % terminate after relative objective value does not changes much.
opts.tol = 10^-5;   % tolerance. 
opts.maxIter = 1000; % maximum iteration number of optimization.
sparsity = zeros(length(lambda), 1);
log_lam  = log(lambda);
for i = 1: length(lambda)
    [W funcVal] = Least_L21(X, Y, lambda(i), opts);
    % set the solution as the next initial point. 
    % this gives better efficiency. 
    opts.init = 1;
    opts.W0 = W;
    sparsity(i) = nnz(sum(W,2 )==0)/d;
end

%%
W_toCsv = [(0:size(W,1)-1)' W];
writematrix(W_toCsv,'W_MTL_Least_L21_MENT14D_STATES.xlsx');
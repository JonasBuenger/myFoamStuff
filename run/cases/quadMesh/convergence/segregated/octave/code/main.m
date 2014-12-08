clear all;
close all;
clc;


%% load simulation data

load('data.mat');

data = cell(4,4);
data(1,:) = {cc1 cv1 s1 T1};
data(2,:) = {cc2 cv2 s2 T2};
data(3,:) = {cc3 cv3 s3 T3};
data(4,:) = {cc4 cv4 s4 T4};
data(5,:) = {cc5 cv5 s5 T5};
%data(6,:) = {cc6 cv6 s6 T6};


%% problem data

Theta_in  = 4;
Theta_out = 1;
R_in      = 0.5;
R_out     = 2.0;
alpha     = 1;


%% calculation

num_cases = 5;
dz = 0.1;
localErrors = cell(num_cases,3);
L2_errors   = zeros(num_cases,3);
num_cells   = zeros(num_cases,1);

for i=1:num_cases
    
    cc = data{i,1};
    cv = data{i,2};
    dA = 1/dz * cv(:,1);
    s_sim = data{i,3};
    T_sim = data{i,4};
    num_cells(i,1) = size(cc,1);
    s_ex = zeros(size(cc,1),2);
    T_ex = zeros(size(cc,1),1);

    [ s_ex(:,1), s_ex(:,2) ] = s_exact( cc, Theta_in, Theta_out, R_in, R_out, alpha );
    [ T_ex(:,1) ] = T_exact( cc, Theta_in, Theta_out, R_in, R_out, alpha );
    
    localErrors{i,1} = s_sim(:,1) - s_ex(:,1);
    localErrors{i,2} = s_sim(:,2) - s_ex(:,2);
    localErrors{i,3} = T_sim(:,1) - T_ex(:,1);

    L2_errors(i,1) = ( sum( localErrors{i,1}.^2 .* dA ) ).^0.5; 
    L2_errors(i,2) = ( sum( localErrors{i,2}.^2 .* dA ) ).^0.5; 
    L2_errors(i,3) = ( sum( localErrors{i,3}.^2 .* dA ) ).^0.5; 

end

%% plot

plotErrors(L2_errors, localErrors, data, num_cells);

%% table of convergence

tableOfConvergence(L2_errors, num_cells)


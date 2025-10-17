clear
close all
clc

Path = ['/home/tanabe/code/magni_data'];
addpath(Path);

%load data file
subject_num = 54;
loadfile = sprintf("dataSeed1_%d",subject_num);
load(loadfile,"p_0","x_0");
load('temp.mat') %サンプリングデータのロード

%for obs 
data = data_ssogmm.xs;


%sampling time
Ts =5;
N = size(data_ssogmm.ts,2);
N_step = 1:N;


%%input
u_m = data_ssogmm.meals;
u_i = data_ssogmm.inses;
modes = data_ssogmm.modes;
p = p_0;
[xe,A_d,B_d1,B_d2,C,i_b,Ge] = linear_matrix(p,Ts);

%%Ra
R_a = get_Ra(data_ssogmm.xs(4,:),p,modes,N);

%%%linear 
delta_x=zeros(5,N);
x=zeros(5,N);

delta_x(:,1) = data_ssogmm.xs([1 2 5 6 7],1)-xe;
x(:,1) = data_ssogmm.xs([1 2 5 6 7],1);

for k=2:N

    delta_u_i = u_i(:,k-1)-i_b;
   
    delta_x(:,k) = A_d*delta_x(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1);
    x(:,k) = delta_x(:,k) + xe;
end


%observer
x_ini = [100;0;5000;1500;300];

L = gain(p,Ts);

delta_xob = zeros(5,N);
xob = zeros(5,N);

delta_xob(:,1) = x_ini - xe;
xob(:,1) = x_ini;

for k=2:N

    delta_y = data(1,k-1) - Ge;
    delta_u_i = u_i(:,k-1) - i_b;

    delta_xob(:,k) = A_d*delta_xob(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1) + L*(delta_y - C*delta_xob(:,k-1));
    xob(:,k) = delta_xob(:,k) + xe;
end


ssogmm = data_ssogmm.xs([1 2 5 6 7],:);
for i=1:5
    figure
    plot(ssogmm(i,:))
    hold on
    plot(x(i,:),'.')
    plot(xob(i,:),'.')
    legend('ssogmm','x','xob')
end


error = ssogmm - xob;
for i=1:5
    figure
    plot(error(i,:))
end

%rmse all
Error_all = rmse(xob',ssogmm');

%rmse only 4day
d = floor(size(N_step,2)/4)*3;
Error_4day = rmse(xob(:,d:end)',ssogmm(:,d:end)');




% % A_pole = eig(A_d)';
% if Ts==5
% %[0.9734 0.9878 0.5432 0.8233 0.8233]%5min
% op =[0.861; 0.88; 0.53; 0.81; 0.52];
% elseif Ts==2.5
% %[0.9866 0.9939 0.7370 0.9073 0.9073]%2.5min
% op = [0.881; 0.93; 0.69; 0.88; 0.7];
% elseif Ts==1
% %[0.9946 0.9976 0.8851 0.9619 0.9619]%1min
% op = [0.91; 0.97; 0.87; 0.95; 0.88];
% end


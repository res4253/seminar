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
p = p_0;

%data
xs = data_ssogmm.xs;
y = data_ssogmm.Gs;

%sampling time
Ts =5;
N = size(data_ssogmm.ts,2);
N_step = 1:N;


%%input
u_m = data_ssogmm.meals;
u_i = data_ssogmm.inses;
modes = data_ssogmm.modes;

%Ra
R_a = get_Ra(xs(4,:),p,modes,N);

U = [u_m;u_i;R_a;modes];


%%
Q = eps; %システム雑音
R = eps; %観測雑音
% v = randn(1,N)*sqrtm(Q);
% w = randn(1,N)*sqrtm(R);


%%ssogmm
dic = dic_ssogmmCL(p,Ts,N,U);
x = dic.dicrete(xs,"ZOH");


%kalman
x_ini = [100;0;5000;1500;300];

% A = @(Ge,Xe) ...
%     [1-Ts*(S_g+Xe), -Ts*Ge, 0, 0, 0;
%     0, 1-Ts*p_2, 0, 0, (Ts*p_2*S_I)/(BW*V_I);
%     0, 0, 1-Ts*k_d, 0, 0;
%     0, 0, Ts*k_d, 1-Ts*k_d, 0;
%     0, 0, 0, Ts*k_d, 1-Ts*k_cl
%     ];
A = dic.A;
A_di = dic.A_di;
A_dx = dic.A_dx;
C = [1 0 0 0 0];

xhat = zeros(5,N);

a = 1;
P = a*eye(5);
xhat(:,1) = x_ini;

b = eye(5);
for k = 2:N

     xhatm = dic.dynamics(xhat(:,k-1),k);
     
     A_tmp = [zeros(3,2) A_di];
     A_temp =[A(xhat(1,k-1),xhat(2,k-1)) ; A_tmp];
     % Pm = A(xhat(1,k-1),xhat(2,k-1))*P*A(xhat(1,k-1),xhat(2,k-1)) + b*Q*b';
     Pm = A_temp*P*A_temp + b*Q*b';

     g = Pm*C'/(C*Pm*C' + R);

     xhat(:,k) = xhatm + g*(y(:,k) - C*xhatm);
     P = (eye(5) - g*C)*Pm;
end


ssogmm = xs([1 2 5 6 7],:);
for  i = 1:5
    figure
    plot(x(i,:));
    hold on
    plot(xhat(i,:));
    plot(ssogmm(i,:))
    legend('x','xaht','ssogmm')
end

error = ssogmm - xhat;
for i=1:5
    figure
    plot(error(i,:))
end

%rmse all
Error_all = rmse(xhat',ssogmm');

%rmse only 4day
d = floor(size(N_step,2)/4)*3;
Error_4day = rmse(xhat(:,d:end)',ssogmm(:,d:end)');
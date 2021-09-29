close all;
clear;clc;
%% ����켣���ٶ�,��ֹ����
%% ˫����
load('..\TrajectoryGenerator\Trajectory_DoubleLaneChange.mat');% [X,Y,theta_des,kappa_des]
Vx=72/3.6;%�̶���������
X_up=140;%�������λ��
%% Alt3_y
% load('..\TrajectoryGenerator\Trajectory_Path_Alt3_y.mat');
% Vx=54/3.6;%�̶���������
% X_up=875;%�������λ��
%% Alt3
% load('..\TrajectoryGenerator\Trajectory_Path_Alt3.mat');
% Vx=54/3.6;%�̶���������
% X_up=875;%�������λ��

%% Straight_plus_CC
% load('..\TrajectoryGenerator\Trajectory_Straight_plus_CC.mat');
% Vx=36/3.6;%�̶���������
% X_up=875;%�������λ��

%% CircleR50
% load('..\TrajectoryGenerator\Trajectory_CircleR50.mat');
% Vx=54/3.6;%�̶���������
% X_up=875;%�������λ��


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Vehicle Parameters for 2DOF and 7DOF                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7DOF&2DOFģ�Ͳ������ṹ��car,����ģ���õĲ���
car.m = 1600;%Mt; % car mass, kg
car.g = 9.8;%g; % gravity acceleration
car.iz = 2000;%Jz; % car moment of inertia about cg, kg m^2
car.lf = 1.016;%a;% distance from front tires to cg, m
car.lr = 1.524;%b; % distance from rear tires to cg, m
car.kyf = 92512;% single front tire cornering stiffness,N/rad 
car.kyr = 67728;% single rear tire cornering stiffness,N/rad
% % ·��Ħ������
car.mu = 1;% ·��Ħ������
% psid_bound=0.85*car.mu*car.g./x0.xd;% 
% beta_bound=atan(0.02*car.mu*car.g); 

car.jw = 1.0;%Jw; % rotor and tire moment of inertia, kg m^2
car.lw = 1.5;%c; % track width, m
car.reff = 0.285;%Rw(1); % tire effective radius, m

%% ������ʼ״̬���ṹ��VehIniSt
VehIniSt.xd = Vx;%V0(1); % initial longitudinal velocity
VehIniSt.yd = 0; % initial lateral velocity
VehIniSt.psid = 0; % initial yaw rate
VehIniSt.omega = VehIniSt.xd/car.reff; % initial wheel rotation speed
VehIniSt.x0 = Tra(1,1);%0; % initial x position
VehIniSt.y0 = Tra(1,2);%0; % initial y position
VehIniSt.psi0 = Tra(1,3);%0; % initial yaw angle

%% �����ó���ģ�Ͳ�����Linear 2DOF
m=car.m;
Iz=car.iz;
lf=car.lf;
lr=car.lr;
l=lf+lr;
Kyf=car.kyf;
Kyr=car.kyr;
% % ��������2DOFģ�͡�����ϵͳPlant����״̬�ռ���
% Ar=[-(2*Kyf+2*Kyr)/m/Vx,-Vx-(2*Kyf*lf-2*Kyr*lr)/m/Vx;
%     -(2*Kyf*lf-2*Kyr*lr)/Iz/Vx,-(2*Kyf*lf.^2+2*Kyr*lr.^2)/Iz/Vx];
% Br=[2*Kyf/m;2*Kyf*lf/Iz];
% Cr=[1,0;0,1];
% Dr=[0;0];

%% (1)�̶�ƫ����㷽ʽ��ѡ������������з��棬������
% �ȹ̶�ƫ�����ģʽ
Error_Mode=0;% ������������һ�����ֵͶӰ�����ƫ��

figure(1);plot(Tra(:,1),Tra(:,2),'r');hold on %�滮�켣

tic
Controller_mode=1;% sfun_LQRwithFF_YC_2DOF
sim('VehicleModelwithLateralController_comp.slx');
toc
figure(1);plot(ans.x,ans.y,'go-','MarkerSize',2);xlabel('X[m]');ylabel('Y[m]');hold on
figure(2);plot(ans.tout,ans.Delta,'go-','MarkerSize',2);xlabel('t[s]');ylabel('Delta[rad]');hold on;

pause;
tic
Controller_mode=2;% sfun_MPC_YC_2DOF
sim('VehicleModelwithLateralController_comp.slx');
toc
figure(1);plot(ans.x,ans.y,'b*-','MarkerSize',2);hold on
figure(2);plot(ans.tout,ans.Delta,'b*-','MarkerSize',2);hold on

pause; 
tic
Controller_mode=3;% sfun_MPC_deltaU_YC_2DOF
sim('VehicleModelwithLateralController_comp.slx');
toc
figure(1);plot(ans.x,ans.y,'md-','MarkerSize',2);hold on
figure(2);plot(ans.tout,ans.Delta,'md-','MarkerSize',2);hold on

pause;
tic
Controller_mode=4;% sfun_MPC_deltaU_softcon_YC_2DOF
sim('VehicleModelwithLateralController_comp.slx');
toc
figure(1);plot(ans.x,ans.y,'k^-','MarkerSize',2);hold on
legend('Target','LQR','MPC','MPC\_dDelta','MPC\_dDelta\_softCon');
figure(2);plot(ans.tout,ans.Delta,'k^-','MarkerSize',2);hold on
legend('LQR','MPC','MPC\_dDelta','MPC\_dDelta\_softCon');

%% (2)�̶���������ѡ��ƫ����㷽ʽ�����з��棬������
% % �ȹ̶�������
% Controller_mode=1;% LQRwithFF
% 
% figure(1);plot(Tra(:,1),Tra(:,2),'r','LineWidth',4);hold on
% tic
% Error_Mode=0;% ������������ƫ��
% sim('VehicleModelwithLateralController_comp.slx');
% toc
% figure(1);plot(ans.x,ans.y,'go-','MarkerSize',2);xlabel('X[m]');ylabel('Y[m]');
% figure(2);plot(ans.x,ans.Delta,'go-','MarkerSize',2);xlabel('X[m]');ylabel('Delta[rad]');hold on;
% pause;
% 
% tic
% Error_Mode=1;Ta=0.2;% ���ݵ�ǰ+Taʱ����λ����Ԥ������ƫ��
% sim('VehicleModelwithLateralController_comp.slx');
% toc
% figure(1);plot(ans.x,ans.y,'bo-','MarkerSize',2);hold on;
% legend('Target','Now','LookAhead');
% figure(2);plot(ans.x,ans.Delta,'bo-','MarkerSize',2);
% legend('Now','LookAhead');

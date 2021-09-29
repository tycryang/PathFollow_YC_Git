function [sys,x0,str,ts,simStateCompliance] = sfun_MPC_YC_2DOF(t,x,u,flag,car,VehIniSt)% (t,x,u,flag,m,Iz,lf,lr,Kyf,Kyr,Vx)
% % 将运算放置于sfun内部，提高独立性
% % 放入初始化模块中，只在仿真开始时运算一次，为了保证变量对output可见，定义了全局变量

switch flag
  case 0   
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(car,VehIniSt);
  case 3
    sys=mdlOutputs(t,x,u);
  case {1,2,4,9}
    sys=[];
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

end


function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(car,VehIniSt)

%% 车辆参数与初始速度
global At Bt Ct Qt Rt lb ub u0
m=car.m;
Iz=car.iz;
lf=car.lf;
lr=car.lr;
Kyf=car.kyf;
Kyr=car.kyr;
Vx=VehIniSt.xd;
% % 基于偏差的状态方程
A=[0,1,0,0;
    0,-(2*Kyf+2*Kyr)/m/Vx,(2*Kyf+2*Kyr)/m,(-2*Kyf*lf+2*Kyr*lr)/m/Vx;
    0,0,0,1;
    0,-(2*Kyf*lf-2*Kyr*lr)/Iz/Vx,(2*Kyf*lf-2*Kyr*lr)/Iz,-(2*Kyf*lf^2+2*Kyr*lr^2)/Iz/Vx];
B=[0;2*Kyf/m;0;2*Kyf*lf/Iz];
% C=[0;-(2*Kyf-2*Kyr)/m/Vx-Vx;-(2*Kyf*lf^2+2*Kyr*lr^2)/Iz/Vx]*psid_ref+[0;0;0;-1]*dpsid_ref;
C=[0;-(2*Kyf-2*Kyr)/m/Vx-Vx;0;-(2*Kyf*lf^2+2*Kyr*lr^2)/Iz/Vx];

% Q=[1 0 0 0;
%     0 0 0 0;
%     0 0 0.05 0;
%     0 0 0 0];%状态权重
Q=[0.05 0 0 0;
    0 0 0 0;
    0 0 1 0;
    0 0 0 0];%状态权重
R=1;%输入权重

% 状态方程离散化
Ts=0.05;
Ad=(eye(4)+A*Ts/2)/(eye(4)-A*Ts/2);
Bd=Ts*B;
Cd=Ts*C;
% MPC预测模型
Np=20; %预测步长，控制步长等于预测步长
At=[]; Bt=[]; Ct=[];
temp1=[];temp2=zeros(size(C,1),1);
Qt=[]; Rt=[];
for i=1:Np
    At=[At; Ad^i];
    Bt=[Bt zeros(size(Bt,1), size(Bd,2));
        Ad^(i-1)*Bd temp1];
    temp1=[Ad^(i-1)*Bd temp1];
    Ct=[Ct;temp2+Ad^(i-1)*Cd];
    temp2=temp2+Ad^(i-1)*Cd;
    
    Qt=[Qt zeros(size(Qt,1),size(Q,1));
        zeros(size(Q,1),size(Qt,1)) Q];
    Rt=[Rt zeros(size(Rt,1),size(R,1));
        zeros(size(R,1),size(Rt,1)) R];
end
% 控制量ut的上下限
lb=-pi/4*ones(Np,1);
ub=pi/4*ones(Np,1);
% 控制量ut的初始值
u0=zeros(Np,1);
    
%% 默认架构
sizes = simsizes;%不动
% 修改以下参数
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
% 修改以上参数
sys = simsizes(sizes);%不动

x0  = [];% initialize the initial conditions
str = [];% str is always an empty matrix
ts  = [Ts 0];% initialize the array of sample times
simStateCompliance = 'UnknownSimState';

end


function sys=mdlOutputs(t,x,u)

global At Bt Ct Qt Rt lb ub u0 %matlab全局变量的要求，在使用的函数中也需要声明一次

err=[u(1),u(2),u(3),u(4)]';
psid_ref=u(5);

% 转换后的优化目标函数矩阵
H=Bt'*Qt*Bt + Rt;
F=Bt'*Qt*(At*err+Ct*psid_ref);
options=optimset('Algorithm','interior-point-convex','Display','off');
coder.extrinsic('quadprog');
ut=quadprog(H,F,[],[],[],[],lb,ub,u0,options);

Delta=ut(1,1);
sys = Delta;

end% end mdlOutputs



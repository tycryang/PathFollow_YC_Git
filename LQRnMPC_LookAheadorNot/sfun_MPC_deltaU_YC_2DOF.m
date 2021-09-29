function [sys,x0,str,ts,simStateCompliance] = sfun_MPC_deltaU_YC_2DOF(t,x,u,flag,car,VehIniSt)% (t,x,u,flag,m,Iz,lf,lr,Kyf,Kyr,Vx)
% % 将运算放置于sfun内部，提高独立性
% % 放入初始化模块中，只在仿真开始时运算一次，为了保证变量对output可见，定义了全局变量
% % 将控制量变为方向盘转角增量

switch flag
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(car,VehIniSt);
  case 3
%       sys=mdlOutputs(t,x,u,At,Bt,Ct,Qt,Rt,Nc,A_cons,b_cons_coeff,dDelta_Max);
    sys=mdlOutputs(t,x,u);
  case {1,2,4,9}
    sys=[];
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

end


function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(car,VehIniSt)

%% 参数初始化
global At Bt Ct Qt Rt Nc A_cons b_cons_coeff dDelta_Max
% % 车辆参数与初始速度
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

Q=[0.05 0 0 0;
    0 0 0 0;
    0 0 1 0;
    0 0 0 0];%状态权重
R=1;%输入权重

% % 状态方程离散化
Ts=0.05;
Ad=(eye(4)+A*Ts/2)/(eye(4)-A*Ts/2);
Bd=Ts*B;
Cd=Ts*C;

% % 扩展状态矩阵，将输入有绝对量变为增量
m = 1;% 控制维数
n = 4;% 原状态维数
Ae = [Ad,Bd;zeros(m,n),eye(m)];
Be = [Bd;eye(m)];
Ce = [Cd;zeros(m)];
Qe = [Q,zeros(n,m);zeros(m,n),R];
Re = R;
% % MPC预测模型
Np=20; %预测步长
Nc=5; % 控制步长
At=[]; Bt_c=[]; Ct=[];
temp1=[];temp2=zeros(size(Ce,1),1);
Qt=[]; Rt=[];
for i=1:Np
    At=[At; Ae^i];
    Bt_c=[Bt_c zeros(size(Bt_c,1), size(Be,2));
        Ae^(i-1)*Be temp1];
    temp1=[Ae^(i-1)*Be temp1];
    Ct=[Ct;temp2+Ae^(i-1)*Ce];
    temp2=temp2+Ae^(i-1)*Ce;
    
    Qt=[Qt zeros(size(Qt,1),size(Qe,1));
        zeros(size(Qe,1),size(Qt,1)) Qe];
end
for j=1:Nc
    Rt=[Rt zeros(size(Rt,1),size(Re,1));
        zeros(size(Re,1),size(Rt,1)) Re];
end
Bt=Bt_c(:,1:Nc);

% % 不等式约束,方向盘转角绝对值上下限
A_k=zeros(Nc,Nc);%见falcone论文 P181
for p=1:Nc
    for q=1:Nc
        if q<=p
            A_k(p,q)=1;
        else
            A_k(p,q)=0;
        end
    end
end
A_I=kron(A_k,eye(m));%对应于falcone论文约束处理的矩阵A,求克罗内克积
A_cons=[A_I; -A_I];
Delta_max=pi/6; % max_steer转角限制
Delta_Max=kron(ones(Nc,1),Delta_max);
b_cons_coeff=[Delta_Max;Delta_Max];
% % 限值约束,方向盘转角速度上下限
% % 这是一种基于舒适性的考虑，限制越严格应该会有越好的舒适性，但是可能会牺牲控制偏差，或者说需要更合理的规划路径
dDelta_max = 0.082;% max_steer_rate转角速度限制 0.0082rad = 0.47deg,
dDelta_Max=kron(ones(Nc,1),dDelta_max);


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

global Delta_aim;
Delta_aim=0;
end


function sys=mdlOutputs(t,x,u)
global At Bt Ct Qt Rt Nc A_cons b_cons_coeff dDelta_Max
global Delta_aim;
err=[u(1),u(2),u(3),u(4),Delta_aim]';
psid_ref=u(5);

% 转换后的优化目标函数矩阵
H=Bt'*Qt*Bt + Rt;
F=Bt'*Qt*(At*err+Ct*psid_ref);

% 完成不等于式约束
Delta_t=kron(ones(Nc,1), Delta_aim);
b_cons=b_cons_coeff+[-Delta_t;Delta_t];
% coder.extrinsic('quadprog');
options=optimset('Algorithm','interior-point-convex','Display','off');
dDelta=quadprog(H,F,A_cons,b_cons,[],[],-dDelta_Max,dDelta_Max,zeros(Nc,1),options);


Delta_aim=Delta_aim +dDelta(1,1); % 取第一个值作为输入
sys = Delta_aim;

end% end mdlOutputs



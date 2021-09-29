function [sys,x0,str,ts,simStateCompliance] = sfun_MPC_YC_2DOF(t,x,u,flag,car,VehIniSt)% (t,x,u,flag,m,Iz,lf,lr,Kyf,Kyr,Vx)
% % �����������sfun�ڲ�����߶�����
% % �����ʼ��ģ���У�ֻ�ڷ��濪ʼʱ����һ�Σ�Ϊ�˱�֤������output�ɼ���������ȫ�ֱ���

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

%% �����������ʼ�ٶ�
global At Bt Ct Qt Rt lb ub u0
m=car.m;
Iz=car.iz;
lf=car.lf;
lr=car.lr;
Kyf=car.kyf;
Kyr=car.kyr;
Vx=VehIniSt.xd;
% % ����ƫ���״̬����
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
%     0 0 0 0];%״̬Ȩ��
Q=[0.05 0 0 0;
    0 0 0 0;
    0 0 1 0;
    0 0 0 0];%״̬Ȩ��
R=1;%����Ȩ��

% ״̬������ɢ��
Ts=0.05;
Ad=(eye(4)+A*Ts/2)/(eye(4)-A*Ts/2);
Bd=Ts*B;
Cd=Ts*C;
% MPCԤ��ģ��
Np=20; %Ԥ�ⲽ�������Ʋ�������Ԥ�ⲽ��
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
% ������ut��������
lb=-pi/4*ones(Np,1);
ub=pi/4*ones(Np,1);
% ������ut�ĳ�ʼֵ
u0=zeros(Np,1);
    
%% Ĭ�ϼܹ�
sizes = simsizes;%����
% �޸����²���
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
% �޸����ϲ���
sys = simsizes(sizes);%����

x0  = [];% initialize the initial conditions
str = [];% str is always an empty matrix
ts  = [Ts 0];% initialize the array of sample times
simStateCompliance = 'UnknownSimState';

end


function sys=mdlOutputs(t,x,u)

global At Bt Ct Qt Rt lb ub u0 %matlabȫ�ֱ�����Ҫ����ʹ�õĺ�����Ҳ��Ҫ����һ��

err=[u(1),u(2),u(3),u(4)]';
psid_ref=u(5);

% ת������Ż�Ŀ�꺯������
H=Bt'*Qt*Bt + Rt;
F=Bt'*Qt*(At*err+Ct*psid_ref);
options=optimset('Algorithm','interior-point-convex','Display','off');
coder.extrinsic('quadprog');
ut=quadprog(H,F,[],[],[],[],lb,ub,u0,options);

Delta=ut(1,1);
sys = Delta;

end% end mdlOutputs



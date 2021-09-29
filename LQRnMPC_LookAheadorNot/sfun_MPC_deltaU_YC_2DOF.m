function [sys,x0,str,ts,simStateCompliance] = sfun_MPC_deltaU_YC_2DOF(t,x,u,flag,car,VehIniSt)% (t,x,u,flag,m,Iz,lf,lr,Kyf,Kyr,Vx)
% % �����������sfun�ڲ�����߶�����
% % �����ʼ��ģ���У�ֻ�ڷ��濪ʼʱ����һ�Σ�Ϊ�˱�֤������output�ɼ���������ȫ�ֱ���
% % ����������Ϊ������ת������

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

%% ������ʼ��
global At Bt Ct Qt Rt Nc A_cons b_cons_coeff dDelta_Max
% % �����������ʼ�ٶ�
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

Q=[0.05 0 0 0;
    0 0 0 0;
    0 0 1 0;
    0 0 0 0];%״̬Ȩ��
R=1;%����Ȩ��

% % ״̬������ɢ��
Ts=0.05;
Ad=(eye(4)+A*Ts/2)/(eye(4)-A*Ts/2);
Bd=Ts*B;
Cd=Ts*C;

% % ��չ״̬���󣬽������о�������Ϊ����
m = 1;% ����ά��
n = 4;% ԭ״̬ά��
Ae = [Ad,Bd;zeros(m,n),eye(m)];
Be = [Bd;eye(m)];
Ce = [Cd;zeros(m)];
Qe = [Q,zeros(n,m);zeros(m,n),R];
Re = R;
% % MPCԤ��ģ��
Np=20; %Ԥ�ⲽ��
Nc=5; % ���Ʋ���
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

% % ����ʽԼ��,������ת�Ǿ���ֵ������
A_k=zeros(Nc,Nc);%��falcone���� P181
for p=1:Nc
    for q=1:Nc
        if q<=p
            A_k(p,q)=1;
        else
            A_k(p,q)=0;
        end
    end
end
A_I=kron(A_k,eye(m));%��Ӧ��falcone����Լ������ľ���A,������ڿ˻�
A_cons=[A_I; -A_I];
Delta_max=pi/6; % max_steerת������
Delta_Max=kron(ones(Nc,1),Delta_max);
b_cons_coeff=[Delta_Max;Delta_Max];
% % ��ֵԼ��,������ת���ٶ�������
% % ����һ�ֻ��������ԵĿ��ǣ�����Խ�ϸ�Ӧ�û���Խ�õ������ԣ����ǿ��ܻ���������ƫ�����˵��Ҫ������Ĺ滮·��
dDelta_max = 0.082;% max_steer_rateת���ٶ����� 0.0082rad = 0.47deg,
dDelta_Max=kron(ones(Nc,1),dDelta_max);


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

global Delta_aim;
Delta_aim=0;
end


function sys=mdlOutputs(t,x,u)
global At Bt Ct Qt Rt Nc A_cons b_cons_coeff dDelta_Max
global Delta_aim;
err=[u(1),u(2),u(3),u(4),Delta_aim]';
psid_ref=u(5);

% ת������Ż�Ŀ�꺯������
H=Bt'*Qt*Bt + Rt;
F=Bt'*Qt*(At*err+Ct*psid_ref);

% ��ɲ�����ʽԼ��
Delta_t=kron(ones(Nc,1), Delta_aim);
b_cons=b_cons_coeff+[-Delta_t;Delta_t];
% coder.extrinsic('quadprog');
options=optimset('Algorithm','interior-point-convex','Display','off');
dDelta=quadprog(H,F,A_cons,b_cons,[],[],-dDelta_Max,dDelta_Max,zeros(Nc,1),options);


Delta_aim=Delta_aim +dDelta(1,1); % ȡ��һ��ֵ��Ϊ����
sys = Delta_aim;

end% end mdlOutputs



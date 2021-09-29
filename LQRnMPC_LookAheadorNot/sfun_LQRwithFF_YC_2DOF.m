function [sys,x0,str,ts,simStateCompliance] = sfun_LQRwithFF_YC_2DOF(t,x,u,flag,car,VehIniSt)% (t,x,u,flag,m,Iz,lf,lr,Kyf,Kyr,Vx)
% % �Ѽ��㲿�ַ���LQRģ���У������㷨�Ķ����Ծͻ�����
% % �Ѽ���ģ�����Iniģ���У���ΪIniģ�鲻����������������Խ�Kc��Delta_ff_coeff����Ϊȫ�ֱ�����ʵ������ʱ��������ͷ���
% % ��ǰ������Ҳ����ģ���У���Ϊǰ����Ҫ�õ�kc3����LQR�Ľ�
%% flag���ú���
switch flag
    case 0
        [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(car,VehIniSt);
    case 3
        sys=mdlOutputs(t,x,u);
    case {1,2,4,9} % unused flags
        sys=[];
    otherwise % error handling
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

end

function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(car,VehIniSt)
global Kc Delta_ff_coeff 
%% ϵͳ������ɢ�������LQR���õ�Kc
% % �����������ʼ�ٶ�
m=car.m;
Iz=car.iz;
lf=car.lf;
lr=car.lr;
Kyf=car.kyf;
Kyr=car.kyr;
Vx=VehIniSt.xd;
% % ����ƫ�������״̬����,״̬=[lateral_error;lateral_error_rate;heading_error;heading_error_rate] 
A=[0,1,0,0;
    0,-(2*Kyf+2*Kyr)/m/Vx,(2*Kyf+2*Kyr)/m,(-2*Kyf*lf+2*Kyr*lr)/m/Vx;
    0,0,0,1;
    0,-(2*Kyf*lf-2*Kyr*lr)/Iz/Vx,(2*Kyf*lf-2*Kyr*lr)/Iz,-(2*Kyf*lf^2+2*Kyr*lr^2)/Iz/Vx];
B=[0;2*Kyf/m;0;2*Kyf*lf/Iz];
% Q=[1 0 0 0;
%     0 0 0 0;
%     0 0 0.05 0;
%     0 0 0 0];% ״̬Ȩ��
% Q=[1 0 0 0;
%     0 0 0 0;
%     0 0 0.5 0;
%     0 0 0 0];% ״̬Ȩ�� %DLC
Q=[0.05 0 0 0;
    0 0 0 0;
    0 0 1 0;
    0 0 0 0];% ״̬Ȩ�� %straightline with CC
R=1;% ����Ȩ��
% % ״̬������ɢ��
Ts=0.05;
Ad=(eye(4)+A*Ts/2)/(eye(4)-A*Ts/2);
Bd=Ts*B;
% % ��ɢLQR����ֵ������
P=Q;
diff=Inf;
i=1;
while(i<1000) && (diff>0.001)
    P_temp=Ad'*P*Ad-Ad'*P*Bd*(R+Bd'*P*Bd)^(-1)*Bd'*P*Ad+Q;% �迨�᷽��/Riccati ����
    diff=max(max(P_temp-P));
    P=P_temp;
    i=i+1;
end
Kc=(R+Bd'*P*Bd)^(-1)*Bd'*P*Ad;

% % ------------- matlab���õ�LQR���� ---------------------
% % dlqr����ɢϵͳ��״̬������������ֵ��Kc�ӽ�
% Kc1=dlqr(Ad,Bd,Q,R);
% % ����ϵͳ״̬��������ɢϵͳ��״̬��������ȱ��Ts����Ϣ����׼�Ȳ���
% Kc2=lqrd(A,B,Q,R,T);
% % ����ϵͳֱ��lqr��״̬������������������ϵͳ
% Kc0 = lqr(A,B,Q,R);
% % ------------------------------------------------------

%% ǰ������
L=lf+lr;
k3=Kc(3);
e2_ss_coeff=-lr+lf*m*Vx.^2/2./Kyr/L;
Kv=m/2/L*(lr/Kyf-lf/Kyr);
Delta_ff_coeff=(L+Kv*Vx.^2+k3*e2_ss_coeff)./Vx;

%% simulinkĬ����Ҫ�Ĳ���
sizes = simsizes;%��Ҫ��
% ������Ĳ���
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
% ������Ĳ���
sys = simsizes(sizes);%��Ҫ��

x0  = [];% initialize the initial conditions
str = [];% str is always an empty matrix
ts  = [Ts 0];% initialize the array of sample times
simStateCompliance = 'UnknownSimState';



end

function sys=mdlOutputs(t,x,u)
global Kc Delta_ff_coeff
err=[u(1);u(2);u(3);u(4)];
psid_ref=u(5);

% ����
Delta_fb = -Kc*err;
% ǰ��
Delta_ff=Delta_ff_coeff*psid_ref;

sys = [Delta_fb;Delta_ff];

end


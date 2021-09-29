function [sys,x0,str,ts,simStateCompliance] = sfun_MPC_deltaU_softcon_YC_2DOF(t,x,u,flag,car,VehIniSt)% (t,x,u,flag,m,Iz,lf,lr,Kyf,Kyr,Vx)
% % �����������sfun�ڲ�����߶�����
% % �����ʼ��ģ���У�ֻ�ڷ��濪ʼʱ����һ�Σ�Ϊ�˱�֤������output�ɼ���������ȫ�ֱ���
% % ����������Ϊ������ת������
% % ������Լ���������ι滮�е��ɳ����Ӻϲ���H����
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
global At Bt Ct Qt Rt Nu Nc A_cons b_cons_coeff lb ub Gay
% % �����������ʼ�ٶ�
m=car.m;
g=car.g;
Iz=car.iz;
lf=car.lf;
lr=car.lr;
Kyf=car.kyf;
Kyr=car.kyr;
mu=car.mu;
Vx=VehIniSt.xd;
% Ħ������״̬����
psid_bound = 0.85*mu*g./Vx;
beta_bound = atan(0.02*mu*g);
ay_bound = mu*g;

% % ������Ӧ����
% % ����ģ�ͣ�״̬Ϊ[beta,wr],Kyf��KyrΪ���ࡢȡ���ĸն�
AAr = [-2*(Kyf+Kyr)/m/Vx,-2*(Kyf*lf-Kyr*lr)/m/Vx^2-1;
    -2*(Kyf*lf-Kyr*lr)/Iz,-2*(Kyf*lf^2+Kyr*lr^2)/Iz/Vx];
BBr = [2*Kyf/m;2*Kyf*lf/Iz];

a11 = AAr(1,1);
a12 = AAr(1,2);
a21 = AAr(2,1);
a22 = AAr(2,2);
b1 = BBr(1);
b2 = BBr(2);
Ta = a11+a22;
Da = a11*a22-a12*a21;

Gr = (b1.*a21-b2.*a11)./Da;%��̬��ڽ��ٶ����棬���뻯������ƫ����Ϊ��̥��ƫ�նȻ��С
Gb = (b2.*a12-b1.*a22)./Da;%��̬���Ĳ�ƫ�����棬���뻯������ƫ����Ϊ��̥��ƫ�նȻ��С
Gay = Gr*Vx;%��̬������ٶ����棬���뻯������ƫ����Ϊ��̥��ƫ�նȻ��С


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
Nu = 1;% ����ά��
Nx = 4;% ԭ״̬ά��
Ae = [Ad,Bd;zeros(Nu,Nx),eye(Nu)];
Be = [Bd;eye(Nu)];
Ce = [Cd;zeros(Nu)];
Qe = [Q,zeros(Nx,Nu);zeros(Nu,Nx),R];
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
A_I=kron(A_k,eye(Nu));%��Ӧ��falcone����Լ������ľ���A,������ڿ˻�
A_cons0=[A_I,zeros(Nu*Nc,1); -A_I,zeros(Nu*Nc,1)];
Delta_max=pi/6; % max_steerת������
Delta_Max=kron(ones(Nc,1),Delta_max);
% b_cons=[Delta_Max-Delta_t;Delta_Max+Delta_t]; %Delta_t��Ҫ��output�����
b_cons_coeff0=[Delta_Max;Delta_Max];
% % ������Լ����Ϊ�������ԣ�ϣ��������ٶ�<4m/s^2���������ӳͷ�
Ay_softbound=4;% Ϊ�������ԣ�ϣ��������ٶ�<4m/s^2���������ӳͷ�
Ay_SoftBound=kron(ones(Nc,1),Ay_softbound);
A_cons1=[Gay*eye(Nu*Nc)*A_I,-ones(Nu*Nc,1);-Gay*eye(Nu*Nc)*A_I,-ones(Nu*Nc,1)];% Gay*Delta_t-epsilon; -Gay*Delta_t-epsilon
b_cons_coeff1=[Ay_SoftBound;Ay_SoftBound];
% % ��������ת��Լ���Ͳ�����ٶ���Լ���ϵ�һ��
A_cons = [A_cons0;A_cons1];
b_cons_coeff = [b_cons_coeff0;b_cons_coeff1];


% % ��ֵԼ��,������ת���ٶ�������,�ɳ�����epsilon������
dDelta_bound = 0.082;% max_steer_rateת���ٶ����� 0.0082rad = 0.47deg
dDelta_Bound=kron(ones(Nc,1),dDelta_bound);
lb = [-dDelta_Bound;0];
ub = [dDelta_Bound;mu*g-Ay_softbound];
% ub = [dDelta_Bound;5*mu*g-Ay_softbound];

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
global At Bt Ct Qt Rt Nu Nc A_cons b_cons_coeff lb ub Gay
global Delta_aim;
err=[u(1),u(2),u(3),u(4),Delta_aim]';
psid_ref=u(5);

% ת������Ż�Ŀ�꺯������
Row=1;
H=[Bt'*Qt*Bt+Rt,zeros(Nu*Nc,1);zeros(1,Nu*Nc),Row];%% Row
F=[Bt'*Qt*(At*err+Ct*psid_ref);0];

% ��ɲ�����ʽԼ��
Delta_t=kron(ones(Nc,1), Delta_aim);
% b_cons=[Delta_Max-Delta_t;Delta_Max+Delta_t];
b_cons=b_cons_coeff+[-Delta_t;Delta_t;-Gay*eye(Nu*Nc)*Delta_t;Gay*eye(Nu*Nc)*Delta_t];
% coder.extrinsic('quadprog');
options=optimset('Algorithm','interior-point-convex','Display','off');
dDelta=quadprog(H,F,A_cons,b_cons,[],[],lb,ub,zeros(Nc+1,1),options);% ��������Լ���ı���

Delta_aim=Delta_aim +dDelta(1,1); % ȡ��һ��ֵ��Ϊ����
sys = Delta_aim;

end% end mdlOutputs



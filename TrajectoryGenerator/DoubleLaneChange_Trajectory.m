clear;
clc;
X=0:0.1:150;
N=length(X);
Y=zeros(1,N);

Par1=fsolve(@Par_dbl1,[0,0,0,0]);
Par2=fsolve(@Par_dbl2,[0,0,0,0]);
for i=1:N
Y(i)=polyval(Par1,X(i))*(X(i)>15 && X(i)<=45)...
    +3.5*(X(i)>45 && X(i)<=70)...
    +polyval(Par2,X(i))*(X(i)>70 && X(i)<=95);
end

X=X';Y=Y';

% % 求曲率和方向角的老方法，因为X等间距，所以方法也对，也没问题
% N=length(Y);
% de_1=diff(Y,1);
% de_2=diff(Y,2);
% Yd=de_1./(X(2)-X(1));
% Ydd=de_2./(X(2)-X(1)).^2;
% Yd(N)=Yd(N-1);
% Ydd(N-1)=Ydd(N-2);
% Ydd(N)=Ydd(N-1);
% psi_ref=atan(Yd);
% kappa_ref=Ydd./(1+Yd.^2).^1.5;

% 路径长度
dX=diff(X);
dY=diff(Y);
ds=sqrt(dX.^2+dY.^2);
s=zeros(length(X),1);
for i=1:length(dX)
    s(i+1)=s(i)+ds(i);
end
% 曲率和方向角
kappa_ref=zeros(N,1);
psi_ref=zeros(N,1);
for i=2:N-1
    x=X(i-1:i+1);
    y=Y(i-1:i+1);
    [kappa_ref(i),psi_ref(i)] = PJcurvature(x,y);
end
kappa_ref(1)=kappa_ref(2);kappa_ref(end)=kappa_ref(end-1);
psi_ref(1)=psi_ref(2);psi_ref(end)=psi_ref(end-1);

% 利用路径长度为横坐标画图
figure(2);set(gcf,'Position',[300,300,900,700]);
subplot(2,2,1);plot(s,X,'o');xlabel('s');ylabel('X');
subplot(2,2,2);plot(s,Y,'o');xlabel('s');ylabel('Y');
subplot(2,2,3);plot(s,psi_ref,'o');xlabel('s');ylabel('theta');
subplot(2,2,4);plot(s,kappa_ref,'o');xlabel('s');ylabel('kappa');
sgtitle('RefTra');

Tra=[X,Y,psi_ref,kappa_ref];
save Trajectory_DoubleLaneChange.mat Tra

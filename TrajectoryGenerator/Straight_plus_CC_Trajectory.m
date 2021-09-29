clear
load('Straight_plus_CC.mat');
X=refx;
Y=refy;
N=length(Y);
% psi_ref=refpsi;

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

% psi_ref=reftheta;
figure(1);plot(X,Y);xlabel('X');ylabel('Y');
% 利用路径长度为横坐标画图
figure(2);set(gcf,'Position',[300,300,900,700]);
subplot(2,2,1);plot(s,X);xlabel('s');ylabel('X');
subplot(2,2,2);plot(s,Y);xlabel('s');ylabel('Y');
subplot(2,2,3);plot(s,psi_ref,'o');xlabel('s');ylabel('psi');
subplot(2,2,4);plot(s,kappa_ref,'o');xlabel('s');ylabel('kappa');
sgtitle('RefTra');

Tra=[X,Y,psi_ref,kappa_ref];
save Trajectory_Straight_plus_CC.mat Tra
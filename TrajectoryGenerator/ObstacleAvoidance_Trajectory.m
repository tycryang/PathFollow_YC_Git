clear;close all;clc;
x=0:0.1:80;
N=length(x);
y=zeros(1,N);
Par1=fsolve(@Par_OA1,[0,0,0,0]);
Par2=fsolve(@Par_OA2,[0,0,0,0]);
for i=1:N
y(1,i)=polyval(Par1,x(i))*(x(i)>12 && x(i)<=25.5)...
    +1*(x(i)>25.5 && x(i)<=36.5)...
    +polyval(Par2,x(i))*(x(i)>36.5 && x(i)<=49);
end
X=x';Y=y';

N=length(Y);
de_1=diff(Y,1);
de_2=diff(Y,2);
Yd=de_1./(X(2)-X(1));
Ydd=de_2./(X(2)-X(1)).^2;
Yd(N)=Yd(N-1);
Ydd(N-1)=Ydd(N-2);
Ydd(N)=Ydd(N-1);
theta_des=atan(Yd);
kappa_des=Ydd./(1+Yd.^2).^1.5;

subplot(3,1,1)
plot(X,Y);
xlabel('X');
ylabel('Y');
subplot(3,1,2)
plot(X,theta_des);
xlabel('X');
ylabel('/theta');
subplot(3,1,3)
plot(X,kappa_des);
xlabel('X');
ylabel('/kappa');

Tra=[X,Y,theta_des,kappa_des]; 
save Trajectory_ObstacleAvoidance.mat Tra
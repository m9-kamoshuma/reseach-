%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
%% Morris-Lecar model

clear;

% {

%% system condition
n=2;

%% ----initial condition----  
t_end=200;
dt=0.001;
tspan=0:dt:t_end;
w=1;
% x0 = 2*w*rand(n,1)-w;  %(-w,w)の範囲でランダム
x0=[-0.2;-0.35];
% x0=[0;-0.1];
% x0=zeros(2,1);

%% ----main program---- %%
[t,x] = ode45(@(t,x) f_ML(x,x0,t),tspan,x0);


%% plot
plot(tspan,x(:,1));hold on;
plot(tspan,x(:,2));
legend('x1', 'x2')


%% モデルの定義
function xdot = f_ML(x,x0,t)

a=1.8;b=3;c=2.2;d=5;
rho=0.3;
u0=0.03;

if t>=100
    u=0.08;
else
    u=0.03;
end

phi1=1/(1+exp(2-4*a*x(1)));
phi2=1/(1+exp(2-4*(d*(x(1)+x0(1)))));

x1dot=c*phi1-b*x(1)-x(2)+u0-u;
x2dot=rho*(phi2-x(2));

    xdot =[x1dot;x2dot];
end

%}

%% calcurate lipsiz constant
%{
X=sym('x',[1 2]) 
f = f_ML(X); 
df = diff(f) 
df = abs(df) 
fplot(df) 

%% モデルの定義
function xdot = f_ML(x)

a=1.8;b=3;c=2.2;d=5;
rho=0.3;
u0=0.03;
u=0;
x0=zeros(2);
% x0=-100*ones(1,2);

phi1=1/(1+exp(2-4*a*x(1)));
phi2=1/(1+exp(2-4*(d*(x(1)+x0(1)))));

x1dot=c*phi1-b*x(1)-x(2)+u0-u;
x2dot=rho*(phi2-x(2));

    xdot =[x1dot;x2dot];
end

%}


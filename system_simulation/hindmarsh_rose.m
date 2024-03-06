%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
%% Hindmarsh-Rose model

clear;
%{
%% system condition
n=3;

%% ----initial condition----  
t_end=1000;
dt=0.001;
tspan=0:dt:t_end;
w=1;
x0 = 2*w*rand(n,1)-w;  %(-w,w)の範囲でランダム
% x0=[0;0;0];

%% ----main program---- %%
[t,x] = ode45(@(t,x) f_hind(x),tspan,x0);


%% plot
plot(tspan,x(:,1));hold on;
plot(tspan,x(:,2));
plot(tspan,x(:,3));
legend('x1', 'x2', 'x3')
%}

%% calcurate lipsiz constant
% {
X=sym('x',[1 3]) 
f = f_quad_hind(X); 
df = diff(f) 
df = abs(df) 
fplot(df) 
%}

%% モデルの定義
function xdot = f_hind(x)
a=1;b=3;c=1;d=5;%定数パラメータ
I=3.25;

r=0.005;s=4;%バースト発火の長さに関わるパラメータ
Q=-1.618;%静止膜電位

x1dot=x(2)-a*x(1)^3+b*x(1)^2-x(3)+I;
x2dot=c-d*x(1)^2-x(2);
x3dot=r*(s*(x(1)-Q)-x(3));

    xdot =[x1dot;x2dot;x3dot];
end

function xdot = f_quad_hind(x)
a=1;b=3;c=1;d=5;%定数パラメータ
I=3.25;

r=0.005;s=4;%バースト発火の長さに関わるパラメータ
Q=-1.618;%静止膜電位

x1dot=x(2)-a*x(1)^3+b*x(1)^2-x(3)+I;
x2dot=c-d*x(1)^2-x(2);
x3dot=r*(s*(x(1)-Q)-x(3));

    xdot =[x(1)x1dot;x2dot;x3dot];
end


%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
clear;

%% system condition
n=3;
N=1;

C=eye(n);
%% ----initial condition----  
t_end=100;
dt=0.001;
tspan=0:dt:t_end;
% tspan=linspace(0,t_end);
b=1;

x0 = 2*b*rand(N*n,1)-b;  %(-b,b)の範囲でランダム
% x0=[0.02;0.01;0.2];
y0 = kron(eye(N),C)*x0;

%% ----main program---- %%
[t,x] = ode45(@(t,x) f_Chua(x),tspan,x0);


%% plot

    plot3(x(:,1),x(:,2),x(:,3));


%% モデルの定義
% ----Chua----

function xdot = f_Chua(x)
%　{
    al=9;%9
    be=14;%14
    ga=0.01;
    k=1;
    A = k*[-al   al    0;
           1     -1    1;
           0    -be   -ga];

    %参考にした文献のfを-fで置き換えているので以下のパラメータ符号が逆になる
    a=1.14;
    b=0.714;
    f=[k*al*(b*x(1)+1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));
        0;
        0];
%}
%{
%   para2
    al=9;
    be=14.28;
    ga=0;
    k=100;
    a=-1/7;%m0
    b=2/7;%m1
    A = k*[-al*b   al    0;
            1     -1    1;
            0    -be   -ga];

    f=-[k*al*(1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));
        0;
        0];
%}

%{
%   para3リミットサイクル？
    al=2;
    be=2;
    ga=0.01;
    k=1;
    a=-1.14;%m0
    b=-0.714;%m1
    A = k*[al*b   al    0;
            1     -1    1;
            0    -be   -ga];

    f=-[k*al*(1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));
        0;
        0];
%}

%Chua parameter4
%{
a=-7.4;
b=-4.1;
c=-1;
k=3.6;
A =    [0   1    0;
           0    0    1;
           a    b     c];
if x(1)>=1/k || x(1)<=-1/k 
    f=[0;
        0;
        k*x(1)];
else 
        f=[0;
        0;
        sign(x(1))];
end
%}

    xdot = A*x + f;
end




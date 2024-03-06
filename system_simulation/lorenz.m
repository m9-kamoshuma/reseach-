%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
%% Morris-Lecar model

clear;

%% system condition
n=3;

%% ----initial condition----  
t_end=50;
dt=0.001;
tspan=0:dt:t_end;
w=1;
x0 = 2*w*rand(n,1)-w;  %(-w,w)の範囲でランダム


%% ----main program---- %%
[t,x] = ode45(@(t,x) f_Lorenz(x),tspan,x0);


%% plot
plot(tspan,x(:,1));hold on;
plot(tspan,x(:,2));
plot(tspan,x(:,3));
legend('x1', 'x2', 'x3')

%% モデルの定義
function xdot = f_Lorenz(x)

sigma=10;r=28;b=8/3;

x1dot=sigma*(x(2)-x(1));
x2dot=r*x(1)-x(2)-x(1)*x(3);
x3dot=x(1)*x(2)-b*x(3);

    xdot =[x1dot;x2dot;x3dot];
end




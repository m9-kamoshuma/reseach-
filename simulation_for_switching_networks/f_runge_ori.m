%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
%4次ルンゲクッタ
function x_d = f_runge_ori(f,x,u,dt)
    k1 = f(x,u)*dt;
    k2 = f(x+0.5*k1,u)*dt;
    k3 = f(x+0.5*k2,u)*dt;
    k4 = f(x+k3,u)*dt;
    x_d=x+(k1+2*k2+2*k3+k4)/6;
end
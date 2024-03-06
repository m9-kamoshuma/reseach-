%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
function xdot = f_Chua_ori(x,u)
% {
%モデルパラメータの設定(from a gallery of Chua attractors)
    al=9;
    be=14;
    ga=0.01;
    k=1;
    A = k*[-al   al    0;
           1     -1    1;
           0    -be   -ga];
    a=1.14;
    b=0.714;
    f=[k*al*(b*x(1)+1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));
        0;
        0];

    xdot = A*x + f + u;
%     xdot = A*x + f;

end




















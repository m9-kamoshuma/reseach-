%% Chua parameter from a gallery of Chua attractors
%para1とpara2でA行列の形と非線形の形が違うため，Gammaの式も変わってくる

%para1(from a gallery of Chua attractors)
%{
al=9;
be=14;
ga=0.01;
k=1;

A = k*[-al   al    0;
       1     -1    1;
       0    -be   -ga];
% 

%nonlinear part
a=1.14;
b=0.714;

% f=[k*al*(b*x(1)+1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));
%     0;
%     0];

Gamma=k*al*max(abs(a),abs(b))*diag([1,0,0]);

%}

%% Chua parameter2
%{
al=9;
be=14.28;
ga=0;
k=0.1;
a=-1/7;%m0
b=2/7;%m1
A = k*[-al*b   al    0;
        1     -1    1;
        0    -be   -ga];
Gamma=k*al*abs(a-b)*diag([1,0,0]);
%}

%% Chua parameter3
%{
al=9;
be=14.28;
ga=0;
k=1;
a=-1/7;%m0
b=2/7;%m1
A = k*[-2.5714   9    0;
        1     -1    1;
        0    -14.286   0];
Gamma=3.8571*diag([1,0,0]);
%}

%% Chua parameter4
% {
al=9;
be=14.28;
ga=0;
k=1;
a=-1/7;%m0
b=2/7;%m1
A = k*[-3.2   10    0;
        1     -1    1;
        0    -15   -0.0385];
Gamma=5.9*diag([1,0,0]);
%}
%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
%%
clear all;
addpath(genpath("C:\yalmip\YALMIP-master"))
addpath(genpath("C:\SeDuMi-master\sedumi-master"))
addpath(genpath("C:\SDPT3-4.0\sdpt3-master"))
addpath(genpath("C:\Program Files\mosek"))
yalmip('clear')

eps = 1e-9;

%% system condition
n=3;

%% Chua parameter from a gallery of Chua attractors
%para1とpara2でA行列の形と非線形の形が違うため，Gammaの式も変わってくる

%para1(from a gallery of Chua attractors)
% {
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
    k=1;
    a=-1/7;%m0
    b=2/7;%m1
    A = k*[-al*b   al    0;
            1     -1    1;
            0    -be   -ga];
Gamma=k*al*abs(a-b)*diag([1,0,0]);
%}
%% other system parameters
B=eye(n);
C=eye(n);
BC=B*C;
normBC=norm(BC);
N=3; %sum of systems


%% LMIparameters
%coupling strength
k=1000;

%% LMI settings
epsN=eps*eye(n);
LMI=[];
P=sdpvar(n,n,'sym');

%% LMI
LMI=[LMI, (A-2*k*BC)'*P + P*(A-2*k*BC)<=epsN];


%% solveLMI
sol=optimize(LMI);
disp(sol.problem);

%% check LMI
%%%%% LMIの値代入%%%%%%%%
P=value(P);
EigP=eig(P);
check_P=min(eig(P));
LMIv=value((A-2*k*BC)'*P + P*(A-2*k*BC));
check_LMI=max(eig(value(LMIv)));









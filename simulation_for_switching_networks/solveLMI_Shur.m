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
%{
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

%% Chua parameter4
% {
a=-7.4;
b=-4.1;
c=-1;
A =    [0   1    0;
           0    0    1;
           a    b     c];
A=-100*ones(n);
Gamma=12*3.6*diag([0,0,1]);
%}
%% other system parameters
B=eye(n);
C=eye(n);
BC=B*C;
normBC=norm(BC);

%% Network conditions
N=3; %sum of systems

%本当はここを可能な数字全てでLMIを解く必要がある
l_st=1; 
l_us=0;

%% LMI
%%%%%%%parameter%%%%%%%%%% %%%
% c=(alpha1+lambda_ast)*T;
% gamma=eta*a2*exp(c)/a1*lambda_ast;
% alpha0 > lambda_ast >lambda
% rate=(alpha1+lambda_ast)/(alpha0-lambda_ast);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%loop results%%
%{
k=100;
alpha0=140;
alpha1=160;
eta=68.9007343514504
taumax=1.00000000000000e-05
lambda_ast=139;
rate=299;
%}
%%%%

k=100;
T=0.00001;
alpha0=300;
alpha1=10000;
lambda_ast=250;
mu=sdpvar(1);%LMI変換のため．mu=1/eta

%%%%%%%%LMI settings%%%%%%%%%%%%%%%%%%
Eps = 1e-9;
EpsN=eps*eye(n);
EpsN2 = eps*eye(n*2);
LMI=[];
P=sdpvar(n,n,'sym');
M22=-inv((N-2)*k^2*(BC)*BC.'+1/Eps*eye(n));

%%%%%%%%%%LMI revised with Shur compliment%%%%%%%%
const_st = [(A-k*(l_st+1)*BC).'*P + P*(A-k*(l_st+1)*BC) + Eps*Gamma+alpha0*P     P;
                   P                                                                                                          M22 ];
const_us = [(A-k*l_us*BC).'*P + P*(A-k*l_us*BC) + Eps*Gamma-alpha1*P                  P;
                   P                                                                                                          M22 ];

LMI=[LMI, const_st <= -EpsN2];
LMI=[LMI, const_us<=-EpsN2];
LMI=[LMI,P>= EpsN];
LMI=[LMI,mu*eye(n)-P<=-EpsN];
LMI=[LMI,mu>=eps];

%固有値の比率a2/a1を１に近づけるためのLMI条件
p1=sdpvar(1);
p2=sdpvar(1);
LMI=[LMI,p1*eye(n)<=P];
LMI=[LMI,P<=p2*eye(n)];

%% solve LMI
ops=sdpsettings('verbose',0);
ops.solver='mosek';

% sol=optimize(LMI,[],ops);
% sol=optimize(LMI);
% sol=optimize(LMI,p2-p1-mu,ops);
sol=optimize(LMI,-mu,ops);

disp(sol.problem);

%% check 
%%%%% LMIの値代入%%%%%%%%
P=value(P);
EigP=eig(P);
check_P=min(eig(P));
const_st=value(const_st);
check_st=max(eig(const_st));
const_us=value(const_us);
check_us=max(eig(const_us));
mu=value(mu);
eta=1/value(mu);
a1=min(eig(P));
a2=max(eig(P));
a2a1=a2/a1;
p1=value(p1);
p2=value(p2);
check_p1LMI=min(eig(P-p1*eye(n)));
check_p2LMI=min(eig(p2*eye(n)-P));



rate=(alpha1+lambda_ast)/(alpha0-lambda_ast);
c=(alpha1+lambda_ast)*T; %2023/6/19に証明して2倍しなくてもよい
gamma=sqrt(((N-2)*eta*a2*exp(c))/(a1*lambda_ast));

%% 表示用
tmp1=sprintf('sol= %d',sol.problem);
tmp2=sprintf('gamma= %d',gamma);
tmp3=sprintf('rate= %d',rate);
tmp4=sprintf('T= %d',T);

disp(tmp1);
disp(tmp2);
disp(tmp3);
disp(tmp4);

if sol.problem==0 && gamma<1 /sqrt(2)
    disp('OK');
end

if sol.problem==0  && check_P>0 && check_st<0 && check_us<0 && eta>0 && check_p1LMI>0 && check_p2LMI>0
    disp('LMI OK');
end








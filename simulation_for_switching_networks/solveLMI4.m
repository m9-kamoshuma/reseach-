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

Gamma=(k*al*max(abs(a),abs(b)))^2*diag([1,0,0]);
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
% A=-eps*eye(n);

%% LMIparameters

k=10;
T=0.00001;
alpha0=20;
alpha1=20;
lambda_ast=11;
eta=5;


%{
k=1000;
T=0.00001;
alpha0=2550;
alpha1=1450;
lambda_ast=1450;
eta=700;
Eps = 1e-6;
%}

%% LMI settings
EpsN=eps*eye(n);
Eps = 1e-6;%ここはLMIの結果にほぼ影響しない
EpsN2= eps*eye(n*2);
z=zeros(n,n);
LMI=[];
P=sdpvar(n,n,'sym');
%1,i間の結合がオンのときの最小の結合数
l_st=1;
%1,i間の結合がオフのときの最小の結合数
l_us=0;

%% LMI

Phi0=(A-k*(l_st+1)*BC).'*P + P*(A-k*(l_st+1)*BC) + Eps*Gamma;
Phi1=(A-k*l_us*BC).'*P + P*(A-k*l_us*BC) + Eps*Gamma;

const_st = [Phi0+alpha0*P    k*P*BC ;    
                   k*BC.'*P          -eta*P    ];

const_us = [Phi1-alpha1*P    k*P*BC;   
                   k*BC.'*P          -eta*P];              

LMI=[LMI, const_st <= -EpsN2];
LMI=[LMI, const_us<=-EpsN2];

% Shur complimentを施す際に必要な条件
const_Shur=(k*BC)'*P + P*(k*BC)-eta*P;
LMI=[LMI,const_Shur>=EpsN];

LMI=[LMI,P>= EpsN];


%固有値の比率a2/a1を１に近づけるためのLMI条件
p1=sdpvar(1);
p2=sdpvar(1);
LMI=[LMI,p1*eye(n)-P<=-EpsN];
LMI=[LMI,P-p2*eye(n)<=-EpsN];


%% solveLMI

ops=sdpsettings('verbose',0);
% ops.solver='mosek';
% sol=optimize(LMI,[],ops);

% sol=optimize(LMI);
sol=optimize(LMI,p2-p1,ops);
disp(sol.problem);

%% check LMI
%%%%% LMIの値代入%%%%%%%%
P=value(P);
EigP=eig(P);
check_P=min(eig(P));
const_st=value(const_st);
check_st=max(eig(const_st));
const_us=value(const_us);
check_us=max(eig(const_us));

a1=min(eig(P));
a2=max(eig(P));
a2a1=a2/a1;


rate=(alpha1+lambda_ast)/(alpha0-lambda_ast);
c=(alpha1+lambda_ast)*T; %2023/6/19に証明して2倍しなくてもよい


gamma=sqrt((2*(N-2)*eta*a2*exp(c))/(a1*lambda_ast));

tmp1=sprintf('sol= %d',sol.problem);
tmp2=sprintf('gamma= %d',gamma);
tmp3=sprintf('rate= %d',rate);
tmp4=sprintf('T= %d',T);

disp(tmp1);
disp(tmp2);
disp(tmp3);
disp(tmp4);

if sol.problem==0 && gamma<1
    disp('OK');
end

if sol.problem==0 && gamma<1 && check_st<0 && check_us<0 && check_P>0
    disp('LMI OK');
end

% for save
%{
progfile=pwd;
datenow = datestr(now,'yyyy-mm-dd-HH-MM');
mkdir('lmi_results',datenow)
cd(strcat('lmi_results/',datenow))
savefile = [sprintf('k%d-',k),sprintf('rate%d-',rate),sprintf('ADT%d-',ADT),datenow,'.mat'];
save(savefile)
cd(progfile)
%}






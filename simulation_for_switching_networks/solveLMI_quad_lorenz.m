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

%% lorenz model

sigma=10;r=28;b=8/3;
% x1dot=sigma*(x(2)-x(1));
% x2dot=r*x(1)-x(2)-x(1)*x(3);
% x3dot=x(1)*x(2)-b*x(3);
A=zeros(n);

%% QUAD condition 
% (Boundedness and synchronization of y-coupled Lorenz systems with or without controllersの(19))
EPS=0.000001;
Br=50;
d=2*b-EPS+(Br-2*(b-EPS)*(sigma-EPS) )^2/(2*Br*(b-EPS))-1;
beta=Br/(2*sigma*(b-EPS))+EPS/sigma;
D=diag([0 d 0]);
P=diag([beta 1 1]);



%% other system parameters
B=eye(n);
C=eye(n);
BC=B*C;
normBC=norm(BC);
N=3; %sum of systems
T=0.00001;%max off dwell time


%% LMI settings
EpsN=eps*eye(n);
Eps = 1e-6;%ここはLMIの結果にほぼ影響しない
EpsN2= eps*eye(n*2);
EpsN3= eps*eye(n*3);
EpsN4= eps*eye(n*4);
EpsN5= eps*eye(n*5);
z=zeros(n,n);
LMI=[];

%1,i間の結合がオンのときの最小の結合数
l_st=1;
%1,i間の結合がオフのときの最小の結合数
l_us=0;

%% LMIparameters
% k=1000;
% T=0.00001;
% alpha0=2550;
% alpha1=1450;
% lambda_ast=1450;
% eta=700;

lambda_ast=100;
k=sdpvar(1);
alpha0=sdpvar(1);
alpha1=sdpvar(1);
eta=sdpvar(1);
LMI=[LMI,k>=eps,alpha0>=lambda_ast,alpha1>=eps,eta>=eps];

%% LMI

Phi0=(A-k*(l_st+1)*BC+D).'*P + P*(A-k*(l_st+1)*BC+D) -2*EPS*eye(n);
Phi1=(A-k*l_us*BC+D).'*P + P*(A-k*l_us*BC+D) -2*EPS*eye(n);


%N=3   
% {
const_st = [Phi0+alpha0*P    k*P*BC   ;
                   k*BC.'*P          -eta*P       ];

const_us = [Phi1-alpha1*P    k*P*BC    ;
                   k*BC.'*P          -eta*P       ];

LMI=[LMI, const_st <= -EpsN2];
LMI=[LMI, const_us<=-EpsN2];
%}

%N=4
%{
const_st = [Phi0+alpha0*P    k*BC*P                k*BC*P           P;
               k*BC.'*P               -eta*P                       z                z; 
               k*BC.'*P                  z                             -eta*P       z;
               P                             z                               z            -Eps*eye(n)];

const_us = [Phi1-alpha1*P    k*BC*P        k*BC*P           P;
               k*BC.'*P          -eta*P                            z                     z; 
               k*BC.'*P               z                            -eta*P               z;
               P                            z                              z                 -Eps*eye(n)];

LMI=[LMI, const_st <= -EpsN4];
LMI=[LMI, const_us<=-EpsN4];
%}

%N=5
%{
const_st = [Phi0+alpha0*P    k*BC*P      k*BC*P     k*BC*P           P;
               k*BC.'*P          -eta*P                   z                     z                     z; 
               k*BC.'*P               z                  -eta*P                 z                    z; 
               k*BC.'*P               z                     z                   -eta*P               z;
               P                            z                     z                     z                 -Eps*eye(n)];

const_us = [Phi1-alpha1*P    k*BC*P      k*BC*P     k*BC*P           P;
               k*BC.'*P          -eta*P                   z                     z                     z; 
               k*BC.'*P               z                  -eta*P                 z                    z; 
               k*BC.'*P               z                     z                   -eta*P               z;
               P                            z                     z                     z                 -Eps*eye(n)];
LMI=[LMI, const_st <= -EpsN5];
LMI=[LMI, const_us<=-EpsN5];
%}



%% solveLMI
ops=sdpsettings('verbose',0);
sol=optimize(LMI,-(alpha0-lambda_ast),ops);

disp(sol.problem);

%% check LMI
%%%%% LMIの値代入%%%%%%%%
k=value(k);
alpha0=value(alpha0);
alpha1=value(alpha1);
eta=value(eta);


const_st=value(const_st);
check_st=max(eig(const_st));
const_us=value(const_us);
check_us=max(eig(const_us));

a1=min(eig(P));
a2=max(eig(P));
a2a1=a2/a1;

rate=(alpha1+lambda_ast)/(alpha0-lambda_ast);
c=(alpha1+lambda_ast)*T; 
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

if sol.problem==0 && gamma<1 && check_st<0 && check_us<0 
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






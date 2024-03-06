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
% {
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


%% LMI
%%%%%%%parameter%%%%%%%%%% 
% c=(alpha1+lambda_ast)*T;
% gamma=eta*a2*exp(c)/a1*lambda_ast;
% alpha0 > lambda_ast >lambda
% rate=(alpha1+lambda_ast)/(alpha0-lambda_ast);

%N=3で可解（新しいgammaの条件は満たしていない．．．）
% k=10000;
% T=0.000001;
% alpha0=19000;
% alpha1=21500;
% lambda_ast=18000;
% eta=10000;

%N=3で可解（新しいgammaの条件は満たしていない．．．）
% k=10000;
% T=0.000001;
% alpha0=19800;
% alpha1=20100;
% lambda_ast=19700;
% eta=10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%
k=50;
T=0.000001;
alpha0=150;
alpha1=50;
lambda_ast=140;
eta=80;

%% LMI settings
EpsN=eps*eye(n);
Eps = 1e-9;
EpsN2= eps*eye(n*2);
EpsN3= eps*eye(n*3);
EpsN6= eps*eye(n*6);
z=zeros(n,n);
LMI=[];
P=sdpvar(n,n,'sym');


%% LMI
%行列をブロック対角にして低次元化する操作をやめて，そのままのLMIを解いた場合

%d_st:(1,k)間結合がオンのときのグラフラプラシアンの次数（結合数）１～N-2
%d_us:(1,k)間結合がオフのときのグラフラプラシアンの次数（結合数）0~N-2
% {
for l_st=1:N-2
    for l_us=0:N-2
        Phi0=(A-k*(l_st+1)*BC).'*P + P*(A-k*(l_st+1)*BC) + Eps*Gamma;
        Phi1=(A-k*l_us*BC).'*P + P*(A-k*l_us*BC) + Eps*Gamma;
        for i2=-1:1
            for i3=-1:1
%                 for i4=-1:1
%                         const_st = [Phi0+alpha0*P    k*i2*BC.'*P      k*i3*BC.'*P     k*i4*BC.'*P           P;
%                                            k*i2*P*BC          -eta*P                   z                     z                     z; 
%                                            k*i3*P*BC               z                  -eta*P                 z                    z; 
%                                            k*i4*P*BC               z                     z                   -eta*P               z;
%                                            P                            z                     z                     z                 -Eps*eye(n)];
%         
%                       const_us = [Phi1-alpha1*P    k*i2*BC.'*P      k*i3*BC.'*P     k*i4*BC.'*P         P;
%                                            k*i2*P*BC          -eta*P                   z                     z                   z; 
%                                            k*i3*P*BC               z                  -eta*P                 z                  z; 
%                                            k*i4*P*BC               z                     z                   -eta*P             z;
%                                            P                            z                     z                     z                -Eps*eye(n)];
%                         LMI=[LMI, const_st <= -EpsN5]
%                         LMI=[LMI, const_us<=-EpsN5]

%N=3を試す場合はendを二つコメントアウト，LMI中のEpsN6をEpsN4に変更       
% {
                        const_st = [Phi0+alpha0*P    k*i2*BC.'*P           P;
                                           k*i2*P*BC          -eta*P                 z; 
                                           P                            z                 -Eps*eye(n)];
        
                      const_us = [Phi1-alpha1*P    k*i2*BC.'*P        P;
                                           k*i2*P*BC          -eta*P              z; 
                                           P                            z           -Eps*eye(n)];
        
                        LMI=[LMI, const_st <= -EpsN3];
                        LMI=[LMI, const_us<=-EpsN3];
%}
        
%                     end
%                 end
            end
        end
    end
end
%}

%{
%%%%%%%大きい行列をシュールコンプリメントした場合％％％％％
for l_st = 1:N-1
    for l_us=0:N-2
            Phi0=(A-k*(l_st+1)*BC).'*P + P*(A-k*(l_st+1)*BC) + Eps*Gamma;
            Phi1=(A-k*l_us*BC).'*P + P*(A-k*l_us*BC) + Eps*Gamma;
            const_st=[Phi0 + alpha0*P + 1/eta*k^2*(N-2)*BC.'*P*BC   P;
                            P                                                           -Eps*eye(n) ];
            const_us=[Phi1 - alpha1*P + 1/eta*k^2*(N-2)*BC.'*P*BC   P;
                            P                                                           -Eps*eye(n) ];
            
            LMI=[LMI, const_st <= -EpsN2];
            LMI=[LMI, const_us<=-EpsN2];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

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


gamma=sqrt(((N-2)*eta*a2*exp(c))/(a1*lambda_ast));

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

if sol.problem==0 && gamma<1/sqrt(2) && check_st<0 && check_us<0 && check_P>0
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






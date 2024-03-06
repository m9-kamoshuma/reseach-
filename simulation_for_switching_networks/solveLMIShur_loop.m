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

%Chua parameter from a gallery of Chua attractors
%{
al=9;
be=14;
ga=0.01;
k=0.1;
A = k*[-al   al    0;
       1     -1    1;
       0    -be   -ga];
%nonlinear part
a=1.14;
b=0.714;
% f=[k*al*(b*x(1)+1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));
%     0;
%     0];
Gamma=k*al*max(abs(a),abs(b))*diag([1,0,0]);
%}

% para2
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

B=eye(n);
C=eye(n);
BC=B*C;
normBC=norm(BC);

%% Network conditions
N=3; %sum of systems
l_st=1;
l_us=0;

%% settings of loops
%k,eta,ltaumax,ambda_astは手探りであらかたの検討をつけて固定する

%coupling strength
k=100;

%converge rate
% alpha0_list=100:100:4000;
alpha0_list=10:10:500;


%diverge rate
% alpha1_list=8000:-200:0;
alpha1_list=1000:-10:10;


%bound of dwell time 
taumax=0.00001;

%bound of converge rate
%lambda_astはalpha0に併せてfor文の中で定義

%for save
solNum=1;
parameter_list=[];

%% LMI
ops=sdpsettings('verbose',0);
EpsN=eps*eye(n);
EpsN2 = eps*eye(n*2);
Eps=1e-9;
M22=-inv((N-2)*k^2*(BC)*BC.'+1/Eps*eye(n));

parfor i_alpha0=1:size(alpha0_list,2)
    alpha0=alpha0_list(i_alpha0);
    for i_alpha1=1:size(alpha1_list,2)
        alpha1=alpha1_list(i_alpha1);
        LMI=[];
        P=sdpvar(n,n,'sym');
        mu=sdpvar(1);%LMI変換のため．mu=1/eta
        p1=sdpvar(1);
        p2=sdpvar(1);

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
        LMI=[LMI,p1*eye(n)<=P];
        LMI=[LMI,P<=p2*eye(n)];
        
        sol=optimize(LMI,p2-p1-mu,ops);
        solNum=sol.problem;
        
        if solNum==0
            %check
            P=value(P);
            eigP=eig(P);      
            check_P=min(eig(P));
            const_st=value(const_st);
            check_st=max(eig(const_st));                        
            const_us=value(const_us);
            check_us=max(eig(const_us));
            mu=value(mu);
            p1=value(p1);
            p2=value(p2);
            check_p1LMI=min(eig(P-p1*eye(n)));
            check_p2LMI=min(eig(p2*eye(n)-P));
  
            if alpha0==140 && alpha1==160
                kkk=0;
            end
            
            if check_P>0 && check_st<0 && check_us<0  &&  check_p1LMI>0 && check_p2LMI>0 && mu>0
                a1=min(eig(P));
                a2=max(eig(P)); 
                eta=1/mu;
           
             % rate
                lambda_ast = alpha0-1;
                rate=(alpha1+lambda_ast)/(alpha0-lambda_ast);
                c=(alpha1+lambda_ast)*taumax;
                gamma=sqrt(((N-2)*eta*a2*exp(c))/(a1*lambda_ast));
                if gamma<1/sqrt(2)
                    parameter=[k;alpha0;alpha1;eta;taumax;lambda_ast;rate];
                    parameter_list=[parameter_list parameter];
                end                            
            end
        end
    end
end

[min_rate,min_rate_index]=min(parameter_list(7,:));
min_rate_para=parameter_list(:,min_rate_index);
          
%% for save
% save('parameter_list.mat','parameter_list')






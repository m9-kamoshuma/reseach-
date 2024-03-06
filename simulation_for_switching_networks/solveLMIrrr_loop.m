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
k=1;
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


%% settings of loops
%k,eta,ltaumax,ambda_astは手探りであらかたの検討をつけて固定する

%coupling strength
k=100;

%converge rate
% alpha0_list=100:100:4000;
alpha0_list=200:1:300;

%diverge rate
% alpha1_list=8000:-200:0;
alpha1_list=300:-5:100;


%gain coefficient
% eta_list=100:100:4000;
eta_list=50:1:100;


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
EpsN2=eps*eye(n*2);
EpsN3=eps*eye(n*3);
EpsN5=eps*eye(n*5);
Eps=1e-6;
l_st=1;
l_us=0;

parfor i_alpha0=1:size(alpha0_list,2)
    alpha0=alpha0_list(i_alpha0);
    for i_alpha1=1:size(alpha1_list,2)
        alpha1=alpha1_list(i_alpha1);
        for i_eta =1:size(eta_list,2)
            eta=eta_list(i_eta);
            if eta>=alpha0/2
                continue
            end

            LMI=[];
            P=sdpvar(n,n,'sym');
            p1=sdpvar(1);
            p2=sdpvar(1);
            z=zeros(n,n);
            
        Phi0=(A-k*(l_st+1)*BC).'*P + P*(A-k*(l_st+1)*BC) + Eps*Gamma;
        Phi1=(A-k*l_us*BC).'*P + P*(A-k*l_us*BC) + Eps*Gamma;
        %N=5
        %{
        const_st = [Phi0+alpha0*P    k*BC*P      k*BC*P     k*BC*P           P;
                       k*BC.'*P          -eta*P                   z                     z                     z; 
                       k*BC.'*P               z                  -eta*P                 z                    z; 
                       k*BC.'*P               z                     z                   -eta*P               z;
                       P                            z                     z                     z                 -Eps*eye(n)];
        
        const_us = [Phi1-alpha1*P    k*i2*BC.'*P      k*i3*BC.'*P     k*i4*BC.'*P         P;
                       k*i2*P*BC          -eta*P                   z                     z                   z; 
                       k*i3*P*BC               z                  -eta*P                 z                  z; 
                       k*i4*P*BC               z                     z                   -eta*P             z;
                       P                            z                     z                     z                -Eps*eye(n)];
        LMI=[LMI, const_st <= -EpsN5]
        LMI=[LMI, const_us<=-EpsN5]
        %}
        
        %N=3   
        % {
        const_st = [Phi0+alpha0*P    k*BC*P           P;
                           k*BC.'*P          -eta*P                 z; 
                           P                            z                 -Eps*eye(n)];
        
        const_us = [Phi1-alpha1*P    k*BC.'*P        P;
                           k*BC.'*P          -eta*P              z; 
                           P                            z           -Eps*eye(n)];
        
        LMI=[LMI, const_st <= -EpsN3];
        LMI=[LMI, const_us<=-EpsN3];

    
            LMI=[LMI,P>= EpsN];
            LMI=[LMI,p1*eye(n)-P<=-EpsN];
            LMI=[LMI,P-p2*eye(n)<=-EpsN];
            
            sol=optimize(LMI,p2-p1,ops);
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
                LMIp1v=value(p1*eye(n)-P);
                check_LMIp1v=max(eig(LMIp1v));
                LMIp2v=value(P-p2*eye(n));
                check_LMIp2v=max(eig(LMIp2v));
                
                
                if check_P>0 && check_st<0 && check_us<0  &&  check_LMIp1v<0 && check_LMIp2v<0
                    a1=min(eig(P));
                    a2=max(eig(P)); 
               
                 % rate
                    lambda_ast = alpha0-1;
                    rate=(alpha1+lambda_ast)/(alpha0-lambda_ast);
                    c=(alpha1+lambda_ast)*taumax;
                    gamma=sqrt(((N-2)*eta*a2*exp(c))/(a1*lambda_ast));
                    if gamma<1/sqrt(2)
                        parameter=[k;alpha0;alpha1;eta;taumax;lambda_ast;rate;gamma];
                        parameter_list=[parameter_list parameter];
                    end                            
                end
            end
        end
    end
end

%{
%search minmum rate
[min_rate,min_rate_index]=min(parameter_list(4,:)/parameter_list(2,:));
min_rate_para=parameter_list(:,min_rate_index);

%search maximum tau
[max_tau,max_tau_index]=max(parameter_list(6,:));
max_tau_para=parameter_list(:,max_tau_index);
               
% for save
save('parameter_list20231115.mat','parameter_list')
%}





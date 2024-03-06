%%%%誤差e_i=x_1-1/N*sum(x_i) %%%%
%%%吉田さんのコードからLMIの部分だけ拝借し，
% サンプリングレートを小さくして疑似的に連続時間としている

clear all;
addpath(genpath("C:\yalmip\YALMIP-master"))
% addpath(genpath("C:\sedumi-master\sedumi-master"))
% addpath(genpath("C:\SDPT3-4.0\sdpt3-master"))
addpath(genpath("C:\Program Files\mosek"))
yalmip('clear')

eps = 1e-9;
%% system
n=3; %dimension
%from a gallery of Chua attractors
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

Gamma=(k*al*max(abs(a),abs(b)))^2*diag([1,0,0]);

% other system parameters
B=eye(n);
C=eye(n);
BC=B*C;
normBC=norm(BC);

%%network 
% L=[1 -1 0; -1 2 -1; 0 -1 1];
L=[1 0 -1; 0  1 -1; -1 -1 2];
% L=[2 -1 -1; -1 1 0; -1 0 1];
% L=[2 -1 -1; -1 2 -1; -1 -1 2];

% L=[1 -1 0; -1 1 0; 0 0 0];
% L=[0 0 0; 0 1 -1; 0 -1 1];

l=eig(L);
N=size(L,1);


c=100; %coupling strength
h=eps;%sampling 

%% LMI settings
EpsN=eps*eye(n);
EpsN2= eps*eye(n*2);
EpsN3= eps*eye(n*3);
EpsN4 = eps*eye(n*4);
EpsN6 = eps*eye(n*6);


%% loop settings
parameter_list=[];
Eps_list = [ 1e-2  1e-6];
lam_list=150:1:200;



% 1onoff2の結果
% lam0=387;lam1=13; rate=29.7

% 1onoff2on3の結果
% lam0=163;lam1=13; rate=12.5

%% LMI
for Epsi=1:size(Eps_list,2)
        for lami = 1:size(lam_list,2)
            for ei=2:N
                Eps=Eps_list(Epsi);
                lam=lam_list(lami);
    
                LMI=[];
                P=sdpvar(n,n,'sym');
                Q=sdpvar(n,n,'sym');
             
                A0 = (A-c*l(ei)*BC);
                B0 = (c*l(ei)*BC);
                M11  = P*A0+A0.'*P + h*A0.'*Q*A0 + Eps*Gamma + lam*P;
                M12  = P*B0 + h*A0.'*Q*B0;
                M13  = P + h*A0.'*Q;
                if lam>0
                    M22  = -(exp(-lam*h)/h)*Q + h*B0.'*Q*B0;
                else
                    M22 = -(1/h)*Q + h*B0.'*Q*B0;
                end
                M23  = h*B0.'*Q;
                M33  = -Eps*eye(n) + h*Q;

                M = [M11      M12     M13;
                         M12.'    M22     M23;
                         M13.'    M23.'   M33];
               
       
                LMI=[LMI,M<=-EpsN3];
                LMI=[LMI,P>=EpsN];
                LMI=[LMI,Q>=EpsN];

                %% solveLMI
                ops=sdpsettings('verbose',0);
                sol = optimize(LMI,[],ops);
                
                %% values
                P=value(P);
                Q=value(Q);               
                M=value(M);
                checkP=min(eig(P));
                checkQ=min(eig(Q));
                checkM=max(eig(M));
                
                if checkP>0 && checkQ>0 && checkM<0 && sol.problem==0 
                    disp('OK')
                    parameter=[Eps;lam;l(ei)];
                    parameter_list=[parameter_list parameter];
                    
                else
                    disp('NO')
                    break;
                end
            end

        end
end








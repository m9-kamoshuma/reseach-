clear all;
addpath(genpath("C:\yalmip\YALMIP-master"))
% addpath(genpath("C:\sedumi-master\sedumi-master"))
% addpath(genpath("C:\SDPT3-4.0\sdpt3-master"))
addpath(genpath("C:\Program Files\mosek"))
yalmip('clear')

eps = 1e-9;
%% system
n=3; %dimension
%from a gallery of Chua attractor
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

%% network 
N=3;
H=[ones(N-1,1) -eye(N-1)];
% pseudo-inverse matrix Hp
Hp=[zeros(1,N-1);-eye(N-1)];
Lon=[1 -1 0; -1 2 -1; 0 -1 1];
Loff=[0 0 0; 0 1 -1; 0 -1 1];
% Loff=[0 0 0; 0 0 0; 0 0 0];
HLonHp=H*Lon*Hp;
HLoffHp=H*Loff*Hp;
c=100; %coupling strength

%% LMI settings
EpsN=eps*eye(n);
EpsN2= eps*eye(n*2);
EpsN4 = eps*eye(n*4);
LMI=[];
P0=sdpvar(n,n,'sym');
P1=sdpvar(n,n,'sym');

%% loop settings
parameter_list=[];
Eps_list = [1  1e-2  1e-4 1e-6];
lam0_list=300:1:500;
lam1_list=30:-1:1;



% 1onoff2の結果
% lam0=387;lam1=13; rate=29.7

% 1onoff2on3の結果
% lam0=163;lam1=13; rate=12.5

%% LMI
parfor Epsi=1:size(Eps_list,2)
    for lam0i = 1:size(lam0_list,2)
        for lam1i = 1:size(lam1_list,2)
            Eps=Eps_list(Epsi);
            lam0=lam0_list(lam0i);
            lam1=lam1_list(lam1i);

            LMI=[];
            P0=sdpvar(n,n,'sym');
            P1=sdpvar(n,n,'sym');

            % 1onoff2
            %{
            phi0=(A-2*c*BC).'*P0 + P0*(A-2*c*BC) + Eps*Gamma + lam0*P0;
            phi1=A.'*P1 + P1*A + Eps*Gamma - lam1*P1;
            
            M0 = [phi0       P0;
                   P0 -Eps*eye(n)];
            M1 = [phi1       P1;
                   P1 -Eps*eye(n)];
            
            LMI=[LMI,M0<=-EpsN2];
            LMI=[LMI,M1<=-EpsN2];
            LMI=[LMI,P0>=EpsN];
            LMI=[LMI,P1>=EpsN];
            %}
            
            % 1onoff2on3
            % {
            phi0=kron(eye(N-1),A.'*P0 + P0*A)-c*kron(HLonHp.',BC.'*P0)...
                -c*kron(HLonHp,P0*BC) + Eps*kron(eye(N-1),Gamma) + lam0*kron(eye(N-1),P0);
            phi1=kron(eye(N-1),A.'*P1 + P1*A)-c*kron(HLoffHp.',BC.'*P1)...
                -c*kron(HLoffHp,P1*BC) + Eps*kron(eye(N-1),Gamma) - lam1*kron(eye(N-1),P1);
            
            M0 =  [phi0         blkdiag(P0,P0);
                  blkdiag(P0,P0)   -Eps*eye(n*(N-1))];
            M1 =  [phi1         blkdiag(P1,P1);
                  blkdiag(P1,P1)   -Eps*eye(n*(N-1))];
            
            LMI=[LMI,M0<=-EpsN4];
            LMI=[LMI,M1<=-EpsN4];
            LMI=[LMI,P0>=EpsN];
            LMI=[LMI,P1>=EpsN];
            %}

            %% solveLMI
            ops=sdpsettings('verbose',0);
            sol = optimize(LMI,[],ops);
            
            %% values
            P0=value(P0);
            P1=value(P1);
            M0=value(M0);
            M1=value(M1);
            checkP0=min(eig(P0));
            checkP1=min(eig(P1));
            checkM0=max(eig(M0));
            checkM1=max(eig(M1));
            
            if checkP0>0 && checkP1>0 && checkM0<0 && checkM1<0 && sol.problem==0 
                disp('OK')
                rate=lam0/lam1;
                parameter=[lam0;lam1;rate];
                parameter_list=[parameter_list parameter];
                
            else
                disp('NO')
                break;
            end

        end
    end
end








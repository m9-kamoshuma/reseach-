%%%%誤差e_i=x_1-1/N*sum(x_i) %%%%

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

%% network 
% N＝３の場合
% {
r=8;
% グラフ読み込み
fname='3sys_all.xlsx';%graph laplasian Excel filename
A0=readmatrix(fname);
N=size(A0,1);
Ag = zeros(N,N,r);
L = zeros(N,N,r);
for i = 1:r
    Ag(:,:,i) = A0(:,(i-1)*N+i:i*N+i-1);
    for j=1:N
        d = sum(Ag(j,:,i));
        L(j,j,i) = d;
    end
    L(:,:,i) = L(:,:,i)-Ag(:,:,i);
end
% ３sysのグラフ番号と形状の対応関係
% 1:完全非連結 
% 2~４:１つon
% 5~7:2つon
% 8:完全グラフ
%}

% N=4の場合
%{
N=4;
r=2^(N*(N-1)/2);
L = zeros(N,N,r);
i=0;
for i1=0:1
    eg1=i1;
    for i2=0:1
        eg2=i2;
        for i3=0:1
            eg3=i3;
            for i4=0:1
                eg4=i4;
                for i5=0:1
                    eg5=i5;
                    for i6=0:1
                        eg6=i6;
                        Ag=[0       eg1   eg5    eg4;
                                eg1   0      eg2    eg6;
                                eg5   eg2    0       eg3;
                                eg4    eg6   eg3       0];

                        i=i+1;
                        for j=1:N
                            d = norm(Ag(j,:),1);
                            L(j,j,i) = d;
                        end
                            L(:,:,i) = L(:,:,i)-Ag;
                    end
                end
            end
        end
    end
end
%}

% N=5の場合
%{
% N=5;
r=1;
L = zeros(N,N,r);
L(:,:,1)=[4 -1 -1 -1 -1;
             -1 3 0 -1 -1;
             -1 0 3 -1 -1;
             -1 -1  -1 4 -1;
             -1 -1 -1 -1 4];
%}



c=100; %coupling strength

%% LMI settings
EpsN=eps*eye(n);
EpsN4 = eps*eye(n*4);
EpsN6 = eps*eye(n*6);
EpsN8 = eps*eye(n*8);
LMI=[];
P=sdpvar(n,n,'sym');

%% loop settings
para_list=[];
para_max_list=[];
Eps_list = [1e-2 1e-6];
dt=1;
lam_list=-14:dt:200;
% lam_list=-50:10:10;

%% LMI
for ri=1:r
    for Epsi=1:size(Eps_list,2)
        Eps=Eps_list(Epsi);
        for lami = 1:size(lam_list,2)
            lam=lam_list(lami);
    
            LMI=[];
            P=sdpvar(n,n,'sym');  

    
            phi=kron(eye(N),A.'*P + P*A)-c*kron(L(:,:,ri),BC.'*P)...
                -c*kron(L(:,:,ri),P*BC) + Eps*kron(eye(N),Gamma) + lam*kron(eye(N),P);
            
            M =  [phi         blkdiag(P,P,P);
                  blkdiag(P,P,P)   -Eps*eye(n*(N))];
          
            LMI=[LMI,M<=-EpsN6];
            LMI=[LMI,P>=EpsN];
    
            %% solveLMI
            ops=sdpsettings('verbose',0);
            sol = optimize(LMI,[],ops);
            
            %% values
            P=value(P);
            M=value(M);
            checkP=min(eig(P));
            checkM=max(eig(M));
            
            if checkP>0  && checkM<0 && sol.problem==0 
    %                 disp('OK')
                para=[ri;Eps;lam];
                para_list=[para_list para];
                
%             else
%                 para_max=[ri;Eps;lam-dt];
%                 para_max_list=[para_max_list para_max];
%                 break;
    %                 disp('NO')
            end
        end
    end

end








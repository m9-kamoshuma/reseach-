%-----3 system-自身-サンプリングのみ-相手-サンプリング後遅延-----%
clear;

%% system condition
n=3;
B=eye(n);
C=eye(n);

%% ----Graph------
N=4;%システム数
r=2;%モード数
% fname='2sysonoff.xlsx';%graph laplasian Excel filename
% fname='3sys(mixed).xlsx';%graph laplasian Excel filename
% fname='12-34.xlsx';%graph laplasian Excel filename
% fname='123-456.xlsx';%graph laplasian Excel filename
fname='4ring-12-34.xlsx';%graph laplasian Excel filename


sheet=1;
L0=readmatrix(fname);
node=size(L0,1);
L = zeros(node,node,r);
for i = 1:r
L(:,:,i) = L0(:,(i-1)*node+i:i*node+i-1);
end
L1=L(:,:,1);
L2=L(:,:,2);

%% -----coupling strength & modes rate-----
%coupling strength
k_max = 20; %20
dk = 0.2; %0.1
k_list = (0:dk:k_max);
k_size = size(k_list,2);

%modes rate
a_max=1;
da=0.01;%0.01
a_list = (0:da:a_max);
a_size = size(a_list,2);

%cycle
cy_max=0.1;
dcy=0.1;
cy_list=(0:dcy:cy_max);
cy_size=size(cy_list,2);


%% ----initial condition----  
NUM=10;
t_end=10;
dt=0.001;%桁数はda*dcy
step=t_end/dt;
t=linspace(0,t_end,step);%similation time
t_size=size(t,2);
b=5;%初期値生成のためのパラメータ
H0=[ones(N-1,1) -eye(N-1)];

%for save
error=zeros(NUM,k_size,a_size,cy_size);
error12=zeros(NUM,k_size,a_size,cy_size);
error13=zeros(NUM,k_size,a_size,cy_size);
error23=zeros(NUM,k_size,a_size,cy_size);
error_ave=zeros(k_size,a_size,cy_size);
error12_ave=zeros(k_size,a_size,cy_size);
error13_ave=zeros(k_size,a_size,cy_size);
error23_ave=zeros(k_size,a_size,cy_size);

synctime=zeros(NUM,k_size,a_size,cy_size);

%% ----main program---- %%
tic
for ki=1:k_size %coupling strength
    k = (ki-1)*dk;
    for ai = 1:a_size %modes rate
        a = (ai-1)*da;
        for cy=1:cy_size 
            cycle=(cy-1)*dcy;
            for num=1:NUM
                x0 = 2*b*rand(N*n,1)-b;  %(-b,b)の範囲でランダム
                y0 = kron(eye(N),C)*x0;
                x=ode4(@(t,x) f_Chua(t,x,N,k,L1,L2,cycle,a),t,x0);
                e=kron(H0,eye(n))*x(end,:).';
                error(num,ki,ai,cy)=norm(e);
                error12(num,ki,ai,cy)=norm(e(1:3));
%                 error13(num,ki,ai,cy)=norm(e(4:6));
                error23(num,ki,ai,cy)=norm(e(1:3)-e(4:6));
            end
        end
    end
end     
toc

for ki=1:k_size %coupling strength
    for ai = 1:a_size %modes rate
        for cy=1:cy_size 
            error_ave(ki,ai,cy)=sum(error(:,ki,ai,cy))/NUM;
            error12_ave(ki,ai,cy)=sum(error12(:,ki,ai,cy))/NUM;
%             error13_ave(ki,ai,cy)=sum(error13(:,ki,ai,cy))/NUM;
            error23_ave(ki,ai,cy)=sum(error23(:,ki,ai,cy))/NUM;           
        end
    end
end     





%% ---for save---
progfile=pwd;
datenow = datestr(now,'yyyy-mm-dd-HH-MM');
mkdir('result_2sysonoff',datenow)
cd(strcat('result_2sysonoff/',datenow))
savefile = [sprintf('NUM%d-',NUM),sprintf('k%d-',max(k_list)),sprintf('cycle%d-',max(cy_list)),sprintf('t_end%d-',t_end),datenow,'.mat'];
save(savefile)
cd(progfile)


%% モデルの定義
% ----Chua----
%xdot：3次元

%モデルパラメータの設定

% chua1
%     A = [2*(0.714-1) 2 0;
%          1 -1 1;
%          0 -2 -0.01];
%     
%     f = [2.28*(abs(x(1)+1)-abs(x(1)-1));
%          0;
%          0];

% chua2
%     A = 0.5[2*(0.714-1) 2 0;
%          1 -1 1;
%          0 -2 -0.01];
%     
%     f = [(abs(x(1)+1)-abs(x(1)-1));
%          0;
%          0];



function xdot = f_Chua(t,x,N,c,L1,L2,cycle,rate)
% paramater chua2
    n=3;

    al=9;%9
    be=14;%14
    ga=0.01;
    k=1;
    A = k*[-al   al    0;
           1     -1    1;
           0    -be   -ga];
    a=-1.14;
    b=-0.714;

    B=eye(n);
    C=eye(n);

  
    f1=-[k*al*(b*x(1)+1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));
    0;
    0];

    f2=-[k*al*(b*x(4)+1/2*(a-b)*(abs(x(4)+1)-abs(x(4)-1)));
    0;
    0];

    if N>=3
        f3=-[k*al*(b*x(7)+1/2*(a-b)*(abs(x(7)+1)-abs(x(7)-1)));
        0;
        0];
    end

    if N>=4
        f4=-[k*al*(b*x(10)+1/2*(a-b)*(abs(x(10)+1)-abs(x(10)-1)));
        0;
        0];
    end
    
    if N==2
        f=[f1;f2];
    elseif N==3
        f=[f1;f2;f3];
    elseif N==4
        f=[f1;f2;f3;f4];
    end
    
    y=kron(eye(N),C)*x;
    if mod(t,cycle)<=cycle*rate
        u = -c.*kron(L1,B)*y;
    else
        u = -c.*kron(L2,B)*y;
    end

    xdot = kron(eye(N),A)*x + f + u;
end
%---------chua_end-----------


%%%%%%%%%
%製作者：加茂脩麻
%%%%%%%%%
clear all;
close all;

%% ---Setting--- %%
% simulation time
T = 0.01;   %simulation time
Tdisp=0.003;
dt = 0.0000001;  %連続時間用
order_dt = 7;
c = 1000;      %結合強度

length=round(T/dt);

%% ----Graph------


%%%%%%グラフ構造が既知の場合%%%%
%{
ver=1;%保存用番号
r=4;
fname='3plus2systems.xlsx';%graph laplasian Excel filename
A0=readmatrix(fname);
N=size(A0,1);
A = zeros(N,N,r);
L = zeros(N,N,r);
% for i = 1:r
% L(:,:,i) = L0(:,(i-1)*node+i:i*node+i-1);
% end
for i = 1:r
    A(:,:,i) = A0(:,(i-1)*N+i:i*N+i-1);
    for j=1:N
        d = sum(A(j,:,i));
        L(j,j,i) = d;
    end
    L(:,:,i) = L(:,:,i)-A(:,:,i);
end
L1=L(:,:,1);%1-2
L2=L(:,:,2);%1-3
L3=L(:,:,3);%2-1-3
L4=L(:,:,4);%no edge
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%未知の場合%%%%%%%%%%%%%%%%
% {
ver=0;%保存用番号
N=3;
Cnum_max=2;%システムの結合最大数
Cnum_min=0;%システムの結合最小数

%初期のグラフ構造
Asum=zeros(1,N);
L=zeros(N);
a=zeros(N);
% a=[0 0 0;0 0 1;0 1 0];
Csum=0;%結合の数
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial condition
n = 3;
b = 1;
x0 = 2*b*rand(N*n,1)-b;  %(-r,r)の範囲でランダム
% x0=[9.1;-5.9;0.1; 5.2;-0.1;5.5; -7.3;8.1;-8.5;];
x=x0;
e_0 = norm(x-kron(ones(N,1),x(1:3,1)));

% input-output
B = eye(n); %input
C = eye(n); %output
y = kron(eye(N),C)*x;

%for switching signal
gsub=zeros(1,N-1);
p=zeros(1,N-1);

% for save
x_t=[]; %state
u_t=[]; %input
y_t=[]; %sampled output
gsum_t=[];
g_t=[];
gsub_t=zeros(N-1,length);
% % a_t=zeros(N,N,length);
% a_t=[0 0 0;0 0 1; 0 1 0];
cycle_t=zeros(100,N-1);
change = zeros(1,N-1);
gsum=1;
%for calc
H0=[ones(N-1,1) -eye(N-1)];
   
%% ---Main program--- %%
i=0;step=1;g=4;
tmpT=zeros(1,N-1);
DToffmax=0.00001;%結合がオフの最大滞留時間(条件値)
DTonmax=0.00003;%結合がオフの最大滞留時間(シミュレーション用設定)
DTon=zeros(1,N-1);%シミュレーション用
DToff=zeros(1,N-1);%シミュレーション用
rate=2.64;%滞留時間のon-offの比率の条件値
rate_supmax=DTonmax/DToffmax;
cycle = zeros(1,N-1);
skiprate=1000;
t_whole=0:dt*skiprate:T;ii=0;x_t_whole=x0;% 全体の描画のためのサンプリング変数
split=(1+rate);da=1/split;%モード切り替え用
for t=0:dt:T

    %Tに対してdtが小さすぎるとき保存する値を間引く
    if mod(i,1000)==0%%
        ii=ii+1;
        x_t_whole(:,ii)=x;
    end
    i=i+1;
    
    % input
    u = -c*kron(L,B)*y;%グラフ未知のとき
%     u = -c*kron(L(:,:,g),B)*y;

    
    % for plot
        x_t(:,i)=x;

% %     y_t(:,i)=y;
% %     u_t(:,i)=u;

        g_t(:,i)=g;

% %     gsum_t(:,i)=Csum;%グラフ未知のとき
    for z=1:N-1
        gsub_t(z,i)=gsub(z);
    end
% %     a_t(:,:,i)=a;%グラフ未知のとき

%  モードの変化記録用  
%{
    if gsub(1)==1
        if gsub(2)==1 && a(2,3)==0
            g=3;
        elseif gsub(2)==1 && a(2,3)==1
            g=4;
        elseif gsub(2)==0 && a(2,3)==0
            g=1;
        else 
           g=2;
        end
    else%gsub(1)==0
        if gsub(2)==1 && a(2,3)==0
            g=5;
        elseif gsub(2)==1 && a(2,3)==1
            g=6;
        elseif gsub(2)==0 && a(2,3)==0
            g=8;
        else 
           g=7;
        end
    end
%}

   
% new state
    for j=1:N
        k=(j-1)*n+1:j*n;
        x(k,1)=f_runge_ori(@f_Chua_ori,x(k,1),u(k,1),dt);

        %外乱を加える
%         if t==0.005
%         x(k,1)=f_runge_ori(@f_Chua_ori,x(k,1),u(k,1),dt)+10*(rand-0.5);
%         end

    end
    
    y=kron(eye(N),C)*x;  

    check_g=gsub;


% 各サブシステムの結合がオフのタイミングをきめる定数をここで決める
    for j=1:N-1

%         切り替えのタイミングのとき
        if mod(t-tmpT(1,j),cycle(j))==0
            change(j) = change(j) + 1;
            tmpT(1,j) = t;%切り替えのタイミングを記録
            DToff(j) = round(rand*DToffmax, order_dt);%DToffの時間を決める
            DTon(j)=round((DTonmax-rate*DToff(j))*rand+rate*DToff(j),order_dt);%DTonの時間を決める
            cycle(j) = DToff(j) + DTon(j);
            cycle_t(change(j),j) = cycle(j);
        end
    
%     結合のon-offを切り替えて，変数を格納   
        if mod(t-tmpT(1,j),cycle(j)) <= DToff(j)
            gsub(j)=0;
        else
            gsub(j)=1;
        end
    end


    % {
%     グラフ構造未知の時
    checkmatrix=(check_g==gsub);
    if min(min(checkmatrix))==0
        A=zeros(N);
        for j=1:N-1
            if gsub(j)==1
                A(1,j+1)=1;
                A(j+1,1)=1;
            end
        end
        tmpA=A;
        Asum=zeros(1,N);
        for a=1:N
            for b=1:N
                Asum(a)=Asum(a)+A(a,b);
            end
        end

        %1以外との結合を追加する
        k=1;
        while k<N%下限なし
%         while min(Asum(2:N))<Cnum_min%結合に下限がある場合
%             if k==N
%                 A=tmpA;
%                 k=0;
%             end

%↓スターグラフにする場合はコメントアウト
            ij=[randi(N-1)+1 randi(N-1)+1];%1との結合はgsubできまるから2からNの乱数
            while ij(1)==ij(2) 
                ij=[randi(N-1)+1 randi(N-1)+1];
            end
           
            if  Asum(ij(1))<Cnum_max && Asum(ij(2))<Cnum_max 
                NUM=randi(2)-1;
                A(ij(1),ij(2))=NUM;
                A(ij(2),ij(1))=NUM;
            end
 %↑スターグラフにする場合はコメントアウト

            Asum=zeros(1,N);
            for a=1:N
                for b=1:N
                    Asum(a)=Asum(a)+A(a,b);
                end
            end
            k=k+1;
        end


        D=zeros(N);
        for z=1:N
            D(z,z)=Asum(z);
        end
        L=D-A;  

        Csum=0;
        for z1=1:N
            for z2=z1+1:N
                Csum=Csum-L(z1,z2);
            end
        end  

        a=A;
    end
    %}

  %グラフ既知のとき 
%{
    if gsub(1)==1 && gsub(2)==0
        g=1;
    elseif gsub(1)==0 && gsub(2)==1
        g=2;
    elseif gsub(1)==1 && gsub(2)==1
        g=3;
    else 
        g=4;
    end
%}
%     ee12=kron(H0,eye(n))*x;
%     e12(i)=norm(ee12(1:2*N));% システム12間の誤差
% 
%     ee23=kron([0 1 -1 0 0],eye(n))*x;
%     e23(i)=norm(ee23);% システム23間の誤差

    e(i)=norm(kron(H0,eye(n))*x);%全体の誤差
end


%% ---Figure settings--- %%
clf
FS = 24;%24
FS_r = 20; % regend 20
FS_a = 24; % axis 24
FONT= 'Times New Roman';
pos1=[0 0 30 30];
pos=[0 0 30 20];
pos2=[0 0 30 10];
pos3=[0 0 30 15];

%学会
pos5=[0 0 15 25];%for sub switching signals
pos6=[0 0 20 20];
FS_r=30;
%学会poster用
pos_poster=[0 0 10 6];


% f1=figure;
% f1.PaperType='a4';
% f1.PaperUnits='centimeters';
% f1.PaperPosition=pos1;
f2=figure;
f2.PaperType='a4';
f2.PaperUnits='centimeters';
f2.PaperPosition=pos_poster;
f3=figure;
f3.PaperType='a4';
f3.PaperUnits='centimeters';
f3.PaperPosition=[0 0 25 8];
f4=figure;
f4.PaperType='a4';
f4.PaperUnits='centimeters';
f4.PaperPosition=pos6;
f5=figure;
f5.PaperType='a4';
f5.PaperUnits='centimeters';
f5.PaperPosition=[0 0 25 8];
% f6=figure;
% f6.PaperType='a4';
% f6.PaperUnits='centimeters';
% f6.PaperPosition=pos2;
f=figure;
f.PaperType='a4';
f.PaperUnits='centimeters';
f.PaperPosition=pos6;
f_state=figure;
f_state.PaperType='a4';
f_state.PaperUnits='centimeters';
f_state.PaperPosition=[0 0 25 8];
f_stateT=figure;
f_stateT.PaperType='a4';
f_stateT.PaperUnits='centimeters';
f_stateT.PaperPosition=[0 0 25 8];
f_error2=figure;


%% 3d state figure

%{
% Fig.1 軌道
figure(f1);
if n == 2
    plot(x_t(1,:),x_t(2,:));
elseif n == 3
    plot3(x_t(1,:),x_t(2,:),x_t(3,:));hold on;
    plot3(x_t(4,:),x_t(5,:),x_t(6,:));hold on;
    plot3(x_t(7,:),x_t(8,:),x_t(9,:));hold on;
    plot3(x_t(10,:),x_t(11,:),x_t(12,:));hold on;
    plot3(x_t(13,:),x_t(14,:),x_t(15,:));
%      parameters = struct('axesPosition', [0.6, 0.1, 0.1, 0.1],...
%                     'zoomZone', [0, 0.001; 0, 1],...
%                     'lineDirection', [1, 2; 4, 3]);
%     zp = BaseZoom();
%     zp.plot(parameters)

end
ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;

xlabel('$$x_1$$','FontSize',FS,'Interpreter','Latex')
ylabel('$$x_2$$','FontSize',FS,'Interpreter','Latex')
zlabel('$$x_3$$','FontSize',FS,'Interpreter','Latex')
%}

%% Chua state figure
%chua単体での状態
%{
figure(f2);
t=0:dt:T;
%N=1;
subplot(n,1,1);
plot(t,x_t(1,:),'LineWidth',1.5); hold on;
ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
% xlabel('time','FontSize',FS,'Interpreter','Latex')
ylabel('$$x_{1,1}$$','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
% lgd=legend('$x_{1,1}$','$x_{2,1}$','$x_{3,1}$','system $x_{2,3}$','system $5$','system $6$','Location','eastoutside','Interpreter','Latex','fontname','Times New Roman');


subplot(n,1,2);
plot(t,x_t(2,:),'LineWidth',1.5); hold on;
ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
% xlabel('time','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
ylabel('$$x_{1,2}$$','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
% legend('$x_{1,2}$','$x_{2,2}$','$x_{3,2}$','system $x_{2,3}$','system $5$','system $6$','Location','eastoutside','Interpreter','Latex','fontname','Times New Roman');


subplot(n,1,3);
plot(t,x_t(3,:),'LineWidth',1.5); hold on;
ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
xlabel('time','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
ylabel('$$x_{1,3}$$','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
% legend('$x_{1,3}$','$x_{2,3}$','$x_{3,3}$','system $x_{2,3}$','system $5$','system $6$','Location','eastoutside','Interpreter','Latex','fontname','Times New Roman');
%}

%% state figure
%{
% Fig.2 各状態
figure(f2);
t=0:dt:T;

tile=tiledlayout(N-1,1,'TileSpacing','none');
for j=1:n
    nexttile
    for i=1:N
        plot(t,x_t((i-1)*n+j,:),'LineWidth',1.5);
        hold on;
    end
    ax = gca;
    ax.FontSize = FS_a;
    ax.FontName = FONT;
    % xlabel('time','FontSize',FS,'Interpreter','Latex')
    ylabel(strcat('$$x',sprintf('_{i,%d}(t)',j),'$$'),'FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
    % lgd=legend('$x_{1,1}$','$x_{2,1}$','$x_{3,1}$','system $x_{2,3}$','system $5$','system $6$','FontSize',FS_r,'Location','eastoutside','Interpreter','Latex','fontname','Times New Roman');
%     zp = BaseZoom();
%     zp.plot;
end
%}

%% state figure
% {
figure(f_state);
sn=1;%図示したい状態の番号
t=0:dt:Tdisp;
for i=1:N
    plot(t,x_t((i-1)*n+sn,1:Tdisp/dt+1),'LineWidth',1.5);
    hold on;
end
 ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
xlabel('time','FontSize',FS,'Interpreter','Latex')
ylabel(strcat('$$x',sprintf('_{i,%d}(t)',sn),'$$'),'FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
% ylim([-1 1])
xlim([0 Tdisp]);
%}

%% whole state figure
% {
figure(f_stateT);
sn=1;%図示したい状態の番号
for i=1:N
        plot(t_whole,x_t_whole((i-1)*n+sn,:),'LineWidth',1.5);
        hold on;
end
ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
xlim([0 T]);
xlabel('time','FontSize',FS,'Interpreter','Latex')
ylabel(strcat('$$x',sprintf('_{i,%d}(t)',sn),'$$'),'FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
%}
%}

%% state 1 in the different scales figure
%{
%上に総時間の状態遷移　下に最初の時間だけの状態遷移
figure(f);
sn=1;%図示したい状態の番号
tiledlayout(2,1,'TileSpacing','tight');
nexttile
for i=1:N
        plot(t_whole,x_t_whole((i-1)*n+sn,:),'LineWidth',1.5);
        hold on;
end
ax = nexttile(1);
ax.FontSize = FS_a;
ax.FontName = FONT;
xlim([0 1]);
% xlabel('time','FontSize',FS,'Interpreter','Latex')
% ylabel(strcat('$$x',sprintf('_{i,%d}(t)',sn),'$$'),'FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
nexttile
t=0:dt:Tdisp;
for i=1:N
    plot(t,x_t((i-1)*n+sn,1:Tdisp/dt+1),'LineWidth',1.5);
    hold on;
end
 ax = nexttile(2);
ax.FontSize = FS_a;
ax.FontName = FONT;
xlabel('time','FontSize',FS,'Interpreter','Latex')
ylabel(strcat('$$x',sprintf('_{i,%d}(t)',sn),'$$'),'FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
% ylim([-1 1])
xlim([0 Tdisp]);
%}



%% error figure
Tdisp=0.005;
% % Fig.3 誤差
% {
figure(f3);
t=0:dt:Tdisp;
plot(t,e(1:Tdisp/dt+1),'b','LineWidth',1.5);
hold on;
% plot(t,e2,'LineWidth',1.5);
% plot(t,e3,'LineWidth',1.5);

ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
xlabel('time','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
ylabel(' $$\|e(t)\|$$','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
xlim([0,Tdisp]);
%}
Tdisp=0.003;
%% 各システムごとの誤差
%{
figure(f_error2);
Tdisp=0.5;
t=0:dt:Tdisp;
plot(t,e12(1:Tdisp/dt+1),'r','LineWidth',1.5);hold on;
plot(t,e23(1:Tdisp/dt+1),'b','LineWidth',1.5);
legend("e12","e23")

ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
xlabel('time','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
ylabel(' $$\|e(t)\|$$','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
xlim([0,Tdisp]);
%}

%% Switching figure

%{
% %Fig.4 スイッチング信号
figure(f4);
hold on;
xlim([0,Tdisp]);
ylim([0.5,8.5]);
% t=0:dt:T;
i=0;
for t=0:dt:Tdisp
    i=i+1;
    plot([t,t+dt],[g_t(i),g_t(i)],"b",'LineWidth',1.5);
    if i~=length+1 && g_t(i)~=g_t(i+1)
        plot([t+dt,t+dt],[g_t(i),g_t(i+1)],"b",'LineWidth',1.5);
    end
end
ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
xlim([0,0.1*Tdisp]);
ylim([0.5,8.5]);
xlabel('time','FontSize',FS,'Interpreter','Latex')
ylabel('Network mode $$\mathcal{M}$$','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
%}

%% Sub switching figure
% {
%Fig.5 スイッチング信号(サブシステム別）
figure(f5);
Tdisp=0.003;
t=tiledlayout(N-1,1,'TileSpacing','none');
for i=2:N
    nexttile
    j=0;
    tmp=0;
    tmpj=1;
    for t=0:dt:Tdisp
        j=j+1;
        if j~=length+1 && gsub_t(i-1,j)~=gsub_t(i-1,j+1)
            plot([tmp,t],[gsub_t(i-1,tmpj),gsub_t(i-1,j)],"b",'LineWidth',1.5);hold on;     
            tmp=t;
            tmpj=j+1;
        end

%         plot([t,t+dt],[gsub_t(i-1,j),gsub_t(i-1,j)],"b",'LineWidth',1.5);hold on;
        if j~=length+1 && gsub_t(i-1,j)~=gsub_t(i-1,j+1)
            plot([t,t],[gsub_t(i-1,j),gsub_t(i-1,j+1)],"b",'LineWidth',1.5);hold on;
        end
    end
    ax = gca;
    if i~=N
        ax.XTick = [];
    end
    ax.FontSize = FS_a;
    ax.FontName = FONT;
    xlim([0,0.1*Tdisp]);
    ylim([-0.5,1.5]);
    ax.YTick = [0 1];
    ylh=ylabel(strcat('$$a',sprintf('_{1 %d}',i),'$$'),'FontSize',FS_a,'Interpreter','Latex','fontname','Times New Roman');
    ylh.Position(2) = ylh.Position(2) - abs(ylh.Position(2) * 0.1);
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle')
end
xlabel('time','FontSize',FS_a,'Interpreter','Latex','fontname','Times New Roman')
%}

%% length of dwelltime 
%{
figure(f7);
t=tiledlayout(N-1,1,'TileSpacing','none');
for i=2:N
    nexttile
    stem(cycle_t(:,i-1));
end
%}

%% coupling figure

%{
% %Fig.6 結合数
figure(f6);
hold on;

i=0;
for t=0:dt:T
    i=i+1;
    plot([t,t+dt],[gsum_t(i),gsum_t(i)],"b",'LineWidth',1.5);
    if i~=length+1 && gsum_t(i)~=gsum_t(i+1)
        plot([t+dt,t+dt],[gsum_t(i),gsum_t(i+1)],"b",'LineWidth',1.5);
    end
end

ax = gca;
ax.FontSize = FS_a;
ax.FontName = FONT;
xlim([0,0.001]);
xlabel('time','FontSize',FS,'Interpreter','Latex')
ylabel('Sum of couplings','FontSize',FS,'Interpreter','Latex','fontname','Times New Roman')
%}



%% make video
%{
%graph movie
j=0;
h=figure;
h.Visible = 'off';
for t=0:dt:T
    j=j+1;
    G = graph(a_t(:,:,j));
    hold off;
    plot(G,'Layout','circle')
    Frame(j) = getframe(h);
end

datenow = datestr(now,'yyyy-mm-dd-HH-MM');
mkdir('video_network',datenow)
cd(strcat('video_network/',datenow))
name = [sprintf('ver%d-',ver),sprintf('N%d-',N),sprintf('Cnum%d-',Cnum_max),sprintf('T%d-',T),'.avi'];
v = VideoWriter(name);
v.FrameRate = 10000; % Framerate

% {
open(v);
writeVideo(v,Frame);
close(v);
%}
%}

%% ---for save---
%{
progfile=pwd;
datenow = datestr(now,'yyyy-mm-dd-HH-MM');
mkdir('simulation_data',datenow)
cd(strcat('simulation_data/',datenow))
name = [sprintf('ver%d-',ver),sprintf('N%d-',N),sprintf('c%d-',c),sprintf('cycle_max%d-',cycle_max),'.mat'];
save(name)
cd(progfile)
%}























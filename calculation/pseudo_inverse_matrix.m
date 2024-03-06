% graph laplacian
% N=3
% L=[2 -1 -1; -1 2 -1; -1 -1 2];%N=3完全グラフ
% L=[1 -1 0; -1 1 0; 0 0 0];%1-2 3
% L=[1 0 -1; 0 0 0; -1 0 1];%1-3 2
% L=[0 0 0; 0 1 -1; 0 -1 1];%1 2-3

% N=4
% L=[3  -1 -1 -1;
%     -1 2  0 -1;
%     -1 0 1 0;
%     -1 -1 0 2];
% L=[3  -1 -1 -1;
%     -1 1  0 0;
%     -1 0 1 0;
%     -1 0 0 1];

% N=4の場合
%{
N=4;
r=2^(N*(N-1)/2);
L = zeros(N,N,r);
A = zeros(N,N,r);
edge_list=[];
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
                            A(:,:,i) = Ag;
                            x=[0 0 1 1];
                            y=[1 0 0 1];
                            plot(graph(Ag),'XData',x,'YData',y,'MarkerSize',5)
                            text(0.3,-0.1,mat2str(round(eig(L(:,:,i)),4)));
                            saveas(gcf,strcat('eig(L)_N=4_',num2str(i),'.png'));
                            edge=[eg1;eg2;eg3;eg4;eg5;eg6];
                            edge=[eg1;eg2;eg3;eg4;eg5;eg6;sum(edge)];
                            edge_list=[edge_list edge];
                    end
                end
            end
        end
    end
end
%}

% N=5
% L=[4 -1 -1 -1 -1; -1 4 -1 -1 -1; -1 -1 3 0 -1; -1 -1 0 3 -1; -1 -1 -1 -1 4];%N=5完全グラフ-(3,4)
% L=[5 -1 -1 -1 -1 -1; 
%     -1 5 -1 -1 -1 -1;
%     -1 -1 4 0 -1 -1;
%     -1 -1  0 4 -1 -1;
%     -1 -1 -1 -1 5 -1;
%     -1 -1 -1 -1 -1 5];
% L=[5 -1 -1 -1 -1 -1; 
%     -1  2 -1  0 0 0;
%     -1 -1  3 -1 0 0;
%     -1  0 -1  2 0 0;
%     -1  0  0  0 2 -1;
%     -1  0  0  0 -1 2];
% 
% L=[5 -1 -1 -1 -1 -1; 
%     -1  3 -1  -1 0 0;
%     -1 -1  4 -1 -1 0;
%     -1  -1 -1  3 0 0;
%     -1  0 -1  0 3 -1;
%     -1  0  0  0 -1 2];


% N=6
% L=[6 -1 -1 -1 -1 -1 -1; 
%     -1 6 -1 -1 -1 -1 -1;
%     -1 -1 5 0 -1 -1 -1;
%     -1 -1  0 5 -1 -1 -1;
%     -1 -1 -1 -1 6 -1 -1;
%     -1 -1 -1 -1 -1 6 -1;
%     -1 -1 -1 -1 -1 -1 6];

N=size(L,1);%number of systems
H=[ones(N-1,1) -eye(N-1)];

% pseudo-inverse matrix Hp
Hp=[zeros(1,N-1);-eye(N-1)];

%% calc
% HLHp=H*L*Hp
% eig(L)
% eig(HLHp)

eigL=zeros(N,1,r);
for i=1:r
    eigL(:,:,i)=eig(L(:,:,i));
end

graph(A);


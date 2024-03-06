% 隣接行列と次数行列からなるグラフラプラシアンのブロック三角化変形

N=4;
W=3;
M=1;
r=[2 1 1];
s=[1 3 4];

A=[0 1  1 0;
      1 0  1 0;
      1 1 0 1;
      0 0 1 0];
D=diag([2,2,3,1]);
% A=[0 1  1;
%     1 0  1;
%       0 1 1];
% D=diag([2,2,2]);

H0=zeros(W,N);
for i=1:W
    for j=1:N
        if j==s(i)
            H0(i,j)=1;
        end
    end
end

for i=1:M
    hi=[ones(r(i)-1,1) -eye(r(i)-1,r(i)-1)];
    if i==1
        H1_left=hi;
    else
        H1_left=blkdiag(H1_left,hi);
    end
end

H1=[H1_left zeros((N-W),2)];

H=[H0;H1];

% 三角化
Dtri=H*D/H;
Atri=H*A/H;
Ltri=Dtri-Atri;


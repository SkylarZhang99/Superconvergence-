clear

% Parameters
Num = 4; % Number of intervals
NO = Num+1; % Number of points
h = 1/Num;

% mesh construction
elem(:,1) = [1:Num];
elem(:,2) = [1:Num]+1;
node = 0:1/Num:1;

DoFs = 2*Num+NO;
bc = 2*Num+[1,NO];
free = setdiff([1:DoFs],bc);

G = (1/Num)*[1/5 1/4 1/3; 1/4 1/3 1/2; 1/3 1/2 1];
B = [-2/3 -1/3 0 1; -1/2 -1/2 0 1; 0 0 -1 1];

At = B'*(G\B);
bt = [h/2; h/2; 0; 0]; % f = 1

w=[1;1];
q = [-1/sqrt(3);1/sqrt(3)]; 
q = [node(:,1:Num) ;node(:,1:Num)]' + h/2 + h/2*repmat(q',Num,1); %matrix for the quadrature points

A = sparse(DoFs,DoFs);
b = sparse(DoFs,1);

for i = 1:Num
    Index = [2*i-1,2*i,2*Num+elem(i,:)];
%     A = A + sparse(ones(4,1)*Index,Index'*ones(1,4),At,DoFs,DoFs);
    A = A + sparse(repmat(Index,4,1),repmat(Index',1,4),At,DoFs,DoFs);
    b = b + sparse(Index,ones(1,4),bt,DoFs,1);
end

% ufree = A(free,free)\b(free);
u = zeros(DoFs,1);
u(free) = A(free,free)\b(free);
% u(bc) = [0; 0];

plot(node,u(2*Num+1:end),'*','LineWidth',2,'MarkerSize',10)% ub
hold on
plot([node(elem(:,2))',node(elem(:,1))']',[u(1:2:2*Num),u(2:2:2*Num)]','-','LineWidth',2);% u0
%% 
exactu=@(x) x.*(1-x)/2;
exactu(q);
eta1=(1-q)/2;
eta2=(1+q)/2
T1=eta1.*exactu(q);
T2=eta2.*exactu(q);
temp1=h*T1*w;
temp2=h*T2*w;
C=[h/3 h/6;h/6 h/3];
S=C\[temp1'; temp2'];
% plot(node,node.*(1-node)/2,'--','LineWidth',2)

% This is Seulip's change.

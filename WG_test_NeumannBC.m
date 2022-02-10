clear all
%format short e
% Parameters
Num = 2; % Number of intervals
NO = Num+1; % Number of points
h = 1/Num;

% mesh construction
elem(:,1) = [1:Num];
elem(:,2) = [1:Num]+1;
node = 0:1/Num:1;
right=node(2:end); 
DoFs = 2*Num+NO;  %number of base functions
bc = 2*Num+[1,NO];  %boundary 
free = setdiff([1:DoFs],bc);   
alpha=0; %neumann boundary condition 
beta=0;  %nuemann boundary condition 


G = h*[1/5 1/4 1/3; 1/4 1/3 1/2; 1/3 1/2 1];
B = [-1/3 -2/3 1 0; -1/2 -1/2 1 0; 0 0 1 -1];


At = B'*(G\B); %local stiffness matrix  
%bt = [h/2; h/2; 0; 0]; % f = 1
Mt = [h/3 h/6 0 0;h/6 h/3 0 0;0 0 0 0; 0 0 0 0];
A = sparse(DoFs,DoFs);
M = sparse(DoFs,DoFs); 
b = sparse(DoFs,1);
n = sparse (DoFs,1); 
n(2*Num+1,1)= -alpha;   %neumann boundary condition 
n(end,1) = beta;       %neumann boundary condition
for i = 1:Num
    Index = [2*i-1,2*i,2*Num+elem(i,:)];
%     A = A + sparse(ones(4,1)*Index,Index'*ones(1,4),At,DoFs,DoFs);
    A = A + sparse(repmat(Index,4,1),repmat(Index',1,4),At,DoFs,DoFs);
   % b = b + sparse(Index,ones(1,4),bt,DoFs,1);
 %     A = A + sparse(ones(4,1)*Index,Index'*ones(1,4),At,DoFs,DoFs);
     M = M + sparse(repmat(Index,4,1),repmat(Index',1,4),Mt,DoFs,DoFs);
end

nq = 7;
[q,w] = lgwt(nq,-1,1);
phi = [(1+q)'/2;(1-q)'/2];
q = repmat(node(:,1:Num),nq,1)' + h/2 + h/2*repmat(q',Num,1);


% computing the right hand side 
f=@(x) pi^2*sin(pi*x);  %f(q) col is qua points row is intervals 
% f=@(x) -(6*x-2);
%f=@(x) -(12*x.^2-12*x+2);
F=h/2*phi*(f(q)'.*w); % check h/2
b(1:2*Num)=F(:);
RHS = b-n;  %neumann boundary condition
% ufree = A(free,free)\b(free);
u = zeros(DoFs,1);
u(end)=1; %fixing one coefficient to guarantee linear independency
u(free) = A(free,free)\RHS(free);  %solution with neumann boundary condition
% u(bc) = [0; 0];

plot(node,u(2*Num+1:end),'*','LineWidth',2,'MarkerSize',10)% ub
hold on
plot([node(elem(:,1))',node(elem(:,2))']',[u(2:2:2*Num),u(1:2:2*Num)]','-','LineWidth',2);% u0 %made changes here: changed [u(1:2:2*Num),u(2:2:2*Num)]'to [u(2:2:2*Num),u(1:2:2*Num)]'
%% 
% exactu=@(x) x.*x.*(x-1);
 exactu=@(x) sin(pi*x);
%exactu=@(x) x.^2.*(x-1).^2;
exactu(q);
eta2=(repmat(right',1,nq)-q)/h;
eta1=ones(Num,nq)-eta2;
T1=eta1.*exactu(q);
T2=eta2.*exactu(q);
temp1=h/2*T1*w;  %h needs to be divided by 2 here by the formula for gaussian quadrature 
temp2=h/2*T2*w;
C=[h/3 h/6;h/6 h/3];
S=C\[temp1'; temp2'];
Qu=zeros(DoFs,1);
Qu(1:2*Num)=S(:);
piu=S(:);
u0=u(1:2*Num);
Qu(2*Num+2:DoFs-1)=exactu(node(2:end-1));
errEn=sqrt((u-Qu)'*A*(u-Qu)); %energy norm need to check order of A , test N-2,4,8,16,32
errL2=sqrt((u-Qu)'*M*(u-Qu));
[errEn errL2]
% plot(node,node.*(1-node)/2,'--','LineWidth',2)
%%
%lifting
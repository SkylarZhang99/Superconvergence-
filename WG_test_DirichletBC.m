clear all
%format short e
% Parameters
% Define problem
TestProblem = 2;
switch TestProblem
    case 1
        exactu = @(x)x.*x.*(x-1);
        f = @(x)-(6*x-2);
        alpha = 0;
        beta = 0;
    case 2
        exactu = @(x)sin(pi*x);
        f=@(x)pi^2*sin(pi*x);
        alpha = 0;
        beta = 0;
    case 3 
        exactu = @(x)x.^2.*(x-1).^2;
        f=@(x) -(12*x.^2-12*x+2);
        alpha = 0;
        beta = 0;
    case 4 
        exactu = @(x)cos(pi*x); 
        f=@(x)pi^2*cos(pi*x);
        alpha = 1;
        beta = -1;
    case 5 
        exactu = @(x)x.^3+x+1;
        f = @(x) -6*x;
        alpha = 1;
        beta = 3;
end

Num = 4; % Number of intervals
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
%alpha=1; %dirichlet boundary condition 
%beta=-1;  %dirichlet boundary condition 


G = h*[1/5 1/4 1/3; 1/4 1/3 1/2; 1/3 1/2 1];
B = [-1/3 -2/3 1 0; -1/2 -1/2 1 0; 0 0 1 -1];
SS = G\B; 
LN = h*[9/5 6/4 1 0; 6/4 4/3 1 0; 1 1 1 0;1/4 1/3 1/2 1];
RN = [3 2 1]'.*G; 

At = B'*(G\B); %local stiffness matrix  
%bt = [h/2; h/2; 0; 0]; % f = 1
Mt = [h/3 h/6 0 0;h/6 h/3 0 0;0 0 0 0; 0 0 0 0];
A = sparse(DoFs,DoFs);
M = sparse(DoFs,DoFs); 
b = sparse(DoFs,1);
d = sparse (DoFs,1); 
d(2*Num+1,1)= alpha;   %dirichlet boundary condition 
d(end,1) = beta;       %dirichlet boundary condition
for i = 1:Num
    Index = [2*i-1,2*i,2*Num+elem(i,:)];
%     A = A + sparse(ones(4,1)*Index,Index'*ones(1,4),At,DoFs,DoFs);
    A = A + sparse(repmat(Index,4,1),repmat(Index',1,4),At,DoFs,DoFs);
   % b = b + sparse(Index,ones(1,4),bt,DoFs,1);
 %     A = A + sparse(ones(4,1)*Index,Index'*ones(1,4),At,DoFs,DoFs);
     M = M + sparse(repmat(Index,4,1),repmat(Index',1,4),Mt,DoFs,DoFs);
end

nq = 2;
[q,w] = lgwt(nq,-1,1);
phi = [(1+q)'/2;(1-q)'/2];
q = repmat(node(:,1:Num),nq,1)' + h/2 + h/2*repmat(q',Num,1);


% computing the right hand side 

F=h/2*phi*(f(q)'.*w); % check h/2
b(1:2*Num)=F(:);
RHS = b-A*d;  %dirichlet boundary condition
% ufree = A(free,free)\b(free);
u = zeros(DoFs,1);
u(2*Num+1)=alpha;
u(end)=beta;
u(free) = A(free,free)\RHS(free);  %solution with dirichlet boundary condition
% u(bc) = [0; 0];
% v = sparse(4,1); % 4 = number of basis functions in one interval

plot(node,u(2*Num+1:end),'*','LineWidth',2,'MarkerSize',10)% ub
hold on
plot([node(elem(:,1))',node(elem(:,2))']',[u(2:2:2*Num),u(1:2:2*Num)]','-','LineWidth',2);% u0 %made changes here: changed [u(1:2:2*Num),u(2:2:2*Num)]'to [u(2:2:2*Num),u(1:2:2*Num)]'
%% 

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
Qu(2*Num+1:DoFs)=exactu(node(1:end));
errEn=sqrt((u-Qu)'*A*(u-Qu)); %energy norm need to check order of A , test N-2,4,8,16,32
errL2=sqrt((u-Qu)'*M*(u-Qu));
[errEn errL2];
% plot(node,node.*(1-node)/2,'--','LineWidth',2)
%%
%lifting
for i = 1:Num
    u_local = [u(2*i-1); u(2*i); u(2*Num+i); u(2*Num+i+1)];
r = SS*u_local;
%r1(1) = r(3);
%r1(2) = r(2);
%r1(3) = r(1);
v = RN*r;
v(4) = h*(u(2*i-1)+u(2*i))/2;

uhat = LN\v;
Uhat(:,i) = uhat;
end
%Uhat = Uhat';
uhatfinal = Uhat(:);
%for i = 1:Num
 %   qn = 7;
  %  [qi,wi] = lgwt(qn,node(i),node(i+1)); 
   % qi = (h*qi + node(i)+node(i+1))/2;
    %yi = @(x) (u(i+1)-u(i))*(x-node(i))/h +u(i);
    %yi = sqrt((exactu(qi) -yi(qi))^2); 
    %yi = yi';
    %l2err(i) = (h*yi*wi)/2;
    
%end
 %l2err = sum(l2err); 
for i = 1:Num
    qn = 7;
    [qi,wi] = lgwt(qn,node(i),node(i+1)); 
    qi = (h*qi + node(i)+node(i+1))/2;
    yii = @(x) uhatfinal(4*i)+(x-node(i))*uhatfinal(4*i-1)/h+(x-node(i)).^2*uhatfinal(4*i-2)/(h^2)+(x-node(i)).^3*uhatfinal(4*i-3)/(h^3);
    yiii = sqrt((exactu(qi) -yii(qi)).^2); 
    yiii = yiii';
    l2errh(i) = (h*yiii*wi)/2;
    xi = node(i):1/2^8:node(i+1);
    plot(xi,yii(xi));
    hold on   
end
   l2errh = sum(l2errh)
   %[l2err   l2errh]
%MMt = [1 1/2 1/3 1/4; 1/2 1/3 1/4 1/5; 1/3 1/4 1/5 1/6; 1/4 1/5 1/6 1/7]; 

%MM = kron(eye(Num),MMt)

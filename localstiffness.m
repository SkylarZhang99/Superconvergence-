function At = localstiffness(node)
Num = 4; % Number of intervals
NO = Num+1; % Number of points
h = 1/Num;
node = sparse ((Num+1)^2,2); 
x = 0:h:1;
y = 0:h:1;
[X,Y] = meshgrid(x,y);
node (:,1) = X(:);
node (:,2) = Y(:);
node = full(node);
elem = sparse (2*Num^2,3);

bc = (Num+1)*(1:Num)';
v = (1:Num*(Num+1))'; 
n = setdiff(v,bc); 
elem (1:Num^2,1) = n; 
elem (Num^2+1:end,1) = n;
elem (1:Num^2,2) = n+Num+1; 
elem (1:Num^2,3) = n+Num+2; 
elem (Num^2+1:end,2) = n+Num+2;    
elem (Num^2+1:end,3) = n+1;          
elem = full(elem);

x1 = node(1,1); y1 = node(1,2);
x2 = node(2,1); y2 = node(2,2);
x3 = node(3,1); y3 = node(3,2);

Ck = [ x1, y1, 1; x2, y2, 1; x3, y3, 1]\eye(3);
Dphi1 = Ck(1:2,1)';
Dphi2 = Ck(1:2,2)';
Dphi3 = Ck(1:2,3)';

area = abs(x1*y2+x2*y3+x3*y1-y1*x2-y2*x3-y3*x1)/2;

At = area*[Dphi1*Dphi1', Dphi2*Dphi1', Dphi3*Dphi1';...
    Dphi1*Dphi2', Dphi2*Dphi2', Dphi3*Dphi2';...
    Dphi1*Dphi3', Dphi2*Dphi3', Dphi3*Dphi3'];

end
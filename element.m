clear all
function [node,elem] = ell(N)
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
end
function At = localstiffness(node)

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
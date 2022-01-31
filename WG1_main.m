clear all
clc

IterNum = 6;
I = (1:IterNum)';
Num = 10*2.^(I-1);

errEn = zeros(IterNum,1);
errL2 = zeros(IterNum,1);

for i = 1:IterNum
    [errEn(i),errL2(i)] = PoissonWG1(Num(i));
end

loglog(Num,errEn,'*-','LineWidth',2,'MarkerSize',10)
hold on
loglog(Num,errEn(end-1)/2*(2^3).^(IterNum-I-1),'--','LineWidth',2)
hold off
title('\fontsize{20}Energy Norm Error')
legend('\fontsize{16}Energy Norm Error','\fontsize{16}Order 3 Convergence')

figure
loglog(Num,errL2,'*-','LineWidth',2,'MarkerSize',10)
hold on
loglog(Num,errL2(end-1)/2*(2^4).^(IterNum-I-1),'--','LineWidth',2)
hold off
title('\fontsize{20}L^2-Norm Error')
legend('\fontsize{16}L^2-Norm Error','\fontsize{16}Order 4 Convergence')
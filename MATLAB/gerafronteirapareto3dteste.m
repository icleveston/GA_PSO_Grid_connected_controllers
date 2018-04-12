a = zeros(2,1);
b = [2;1];
c = [3;-.5];
fun = @(u)simple_multiobj2(u,a,b,c);
[x,f,ef] = gamultiobj(fun,2)
figure(1)
scatter3(f(:,1),f(:,2),f(:,3),'k.');

figure(2)
F = scatteredInterpolant(f(:,1),f(:,2),f(:,3),'linear','none');
sgr = min(f(:,1)):.01:max(f(:,1));
ygr = min(f(:,2)):.01:max(f(:,2));
[XX,YY] = meshgrid(sgr,ygr);
ZZ = F(XX,YY);
surf(XX,YY,ZZ,'LineStyle','none')
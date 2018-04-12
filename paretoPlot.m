load('pareto.mat');
f = pareto';
f(:,2) = f(:,2)/1e4;
F = scatteredInterpolant(f(:,1),f(:,2),f(:,3),'linear','none');
sgr = min(f(:,1)):.01:max(f(:,1));
ygr = min(f(:,2)):.01:max(f(:,2));
[XX,YY] = meshgrid(sgr,ygr);
ZZ = F(XX,YY);

scatter3(0.11578, 14711.82162/1e4, 0.99778, 'MarkerEdgeColor','b', 'MarkerFaceColor','b');
hold;
mesh(XX,YY,ZZ, 'EdgeColor', 'b');
xlabel('\beta');
ylabel('\epsilon');
zlabel('\sigma');
function f = simple_multiobj2(x,a,b,c)
x = x(:); a = a(:); b = b(:); c = c(:); % all column vectors
f(1) = sqrt(1+norm(x-a)^2);
f(2) = 0.5*sqrt(1+norm(x-b)^2) + 2;
f(3) = 0.25*sqrt(1+norm(x-c)^2) - 4;


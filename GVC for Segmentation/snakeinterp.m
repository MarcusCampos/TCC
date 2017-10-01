function [xi,yi] = snakeinterp(x,y,dmax,dmin)

x = x(:); y = y(:);  num=0;
N = length(x);
d = abs(x([2:N 1])- x(:)) + abs(y([2:N 1])- y(:));

%   d(i,i+1)<dmin, then either i or i+1 is removed 
IDX = (d<dmin);
idx = find(IDX==0);
x = x(idx);
y = y(idx);
N = length(x);
d = abs(x([2:N 1])- x(:)) + abs(y([2:N 1])- y(:));


%   d(i,i+1)>dmax, then a new point is added between i and i+1
IDX = (d>dmax);
z = snakeindex(IDX);  % interpolate index
p = 1:N+1;
xi = interp1(p,[x;x(1)],z');
yi = interp1(p,[y;y(1)],z');

N = length(xi);
d = abs(xi([2:N 1])- xi(:)) + abs(yi([2:N 1])- yi(:));

while (max(d)>dmax),

    IDX = (d>dmax);
    z = snakeindex(IDX);

    p = 1:N+1;

    xi = interp1(p,[xi;xi(1)],z');
    yi = interp1(p,[yi;yi(1)],z');

    N = length(xi);
    d = abs(xi([2:N 1])- xi(:)) + abs(yi([2:N 1])- yi(:));
end

function [x,y] = snakeinit(delta)
hold on
x = [];
y = [];
n =0;

disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')

but = 1;
while but == 1
      [s, t, but] = ginput(1);
      n = n + 1;
      x(n,1) = s;
      y(n,1) = t;
      plot(x, y, 'r-');
end   

plot([x;x(1,1)],[y;y(1,1)],'r-');
hold off


% linear interpolation
x = [x;x(1,1)];
y = [y;y(1,1)];
t = 1:n+1;
ts = [1:delta:n+1]';
xi = interp1(t,x,ts);
yi = interp1(t,y,ts);
n = length(xi);
x = xi(1:n-1);
y = yi(1:n-1);

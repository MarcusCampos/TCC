function snakedisp(x,y,style)
hold on
x = x(:); y = y(:);
plot([x;x(1,1)],[y;y(1,1)],style);
hold off

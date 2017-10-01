function   [x,y] = reAlign(x,y)
N = length(x);
xbar = mean( x ); ybar = mean( y );
Z = (x - xbar) + (y - ybar)*i; %complex
theta = angle(Z) ;% theta [-pi,+pi]
nn = find(theta < 0.0);
theta(nn) = theta(nn) + 2*pi; % [-pi,0]-->[pi,2*pi], theta [0,+2*pi]

for j = 1:N,
    [mintheta, idx] = min(theta);
    xtemp(j) = x(idx);
    ytemp(j) = y(idx);
    theta(idx) = 2*pi + theta(idx);% add 2*pi after ranking, make theta bigger than 2*pi
end
 x = xtemp';
 y = ytemp';

function [X0,Y0] = alignEndoToEpi(len,X00,Y00)
%
x = X00;
y = Y00;
N = length(X00);

if len == N,
    X0 = X00;
    Y0 = Y00;
   return;
else
    d = sqrt( ( x([2:N 1])- x(:) ).^2.0 + ( y([2:N 1])- y(:) ).^2.0 ) ;
    
    dd = d;
    maxNo = 0;
%  & nbsp;  
    No = len - N;
    for i = 1:No,% 
        [deltaD,indx] = max(dd);  % search the max distance adjacent points
        dd(indx) = 0;
        maxNo = maxNo + length(indx);
        if maxNo >= len-N, break;end
    end
    
    IDX = (d>=deltaD);
    noIDX = sum(IDX);
    if  noIDX < No,
        error('There is no enough points to interpolate......');
    elseif noIDX >No,
        count = 0;
        for j = 1:N,
            if IDX(j) == 1,
                IDX(j)= 0;
%       &nbsp;         
                count = count +1;
            end
            if count == noIDX -No, break;end
        end
    end
    
    i = 1:N;
    flagf(2*i -1) = 0;
    flagf(2*i) = 0;
    flagf(2*i(IDX == 0)) = [];
    
    z = snakeindex(IDX);
    
    p = 1:N+1;
    
    xi = interp1(p,[x;x(1)],z');
    yi = interp1(p,[y;y(1)],z');
    
    %---adjust the spacing-------
    A = CreateCoefMatrix(len, 0,1);
    invAI = inv(A + eye(len));
    xi = invAI * ((xi) );%
    yi = invAI * ((yi) );
% & nbsp;   

    %---adjust end---------------
    
    [X0,Y0] = reAlign(xi,yi);
    if length(X0) ~= len,error('Error in interpolation.......');end
    return;
end
end



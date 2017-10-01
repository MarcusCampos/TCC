function f = NewEdgeMap(Ima, x,y)

[ffx,ffy] = gradient(Ima); % Calculate the gradient of the edge map
ff = sqrt(ffx.*ffx+ffy.*ffy);

xx = round(x);
yy = round(y);
m = length(x);
b = max(yy);
a = min(yy);
for t = a : b
    k = 1;
    for i = 1 : m%scan the edge point         
        if yy(i) == t
            aa(k) = xx(i);
            k = k+1;
        end
    end   
    for i = min(aa(:)) : max(aa(:))
        ff(t,i) = 0;
    end
end 

[fx,fy]=gradient(ff);
f = fx.^2 + fy.^2;
end


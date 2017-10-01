function Y = xconv2(I,G)
[n,m] = size(I);
[n1,m1] = size(G);
FI = fft2(I,n+n1-1,m+m1-1);  % avoid aliasing
FG = fft2(G,n+n1-1,m+m1-1);
FY = FI.*FG;
YT = real(ifft2(FY));
nl = floor(n1/2);
ml = floor(m1/2);
Y = YT(1+nl:n+nl,1+ml:m+ml);

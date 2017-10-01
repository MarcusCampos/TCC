function [fx, fy] = MyGVC(f,rx,ry,h,alpha)
% reference:
%  Hyun Keun Park and Myung Jin Chung, External force of snake: virtual electric field,IEE Electronics Letters,38(24)1500-1502, 2002
%  Dan Yuan, Siwei Lu, Simulated static electric field (SSEF) snake for deformable models,IEEE ICIP, 83-86, 2002
%  Zhanjun Yue, Ardeshir Goshtasby, Laurens V Ackerman, Automatic detection of rib borders in chest radiographs, IEEE trans Med. Imag., 14(3):525-536, 1995
%  Andrei C. Jalba, M.H.F. Wilkinson, and J.B.T.M. Roerdink, CPM: A Deformable Model for Shape Recovery and Segmentation Based on Charged Particles, IEEE TPAMI, 26(10) 1320-1335, 2004
%  Andrew DJ Cross, Edwin R. Hancock, scale space vector fields for symmetry detection, IVC, 17:337-345,1999

% time=clock;
[m,n] = size(f);
fmin  = min(f(:));
fmax  = max(f(:));
f = (f-fmin)/(fmax-fmin); 

[Mx,My] = createMask(rx,ry,h,alpha);
fx = xconv2(f,Mx);
fy = xconv2(f,My);

% fprintf(1,'  in GVC, rx [%d], ry [%d]\n', rx,ry);
% time = etime(clock,time)

function [Mx,My] = createMask(rx,ry,h,alpha)
Rx = floor(rx) - 1;
Ry = floor(ry) - 1;

for i = -Ry:Ry,
    for j = -Rx:Rx,
        if i == 0 & j == 0,
            Mx(i+ Ry+1,j+Rx+1) = 0;
            My(i+ Ry+1,j+Rx+1) = 0;           
            continue;
        end
        Mx(i+ Ry+1,j+Rx+1) = -j/power(sqrt(i*i + j*j + h),alpha);
        My(i+ Ry+1,j+Rx+1) = -i/power(sqrt(i*i + j*j + h),alpha);
        
    end
end

 function Y = xconv2(I,G)
[n,m] = size(I);
[n1,m1] = size(G);
FI = fft2(I,n+n1-1,m+m1-1);  
FG = fft2(G,n+n1-1,m+m1-1);
FY = FI.*FG;
YT = real(ifft2(FY));
nl = floor(n1/2);
ml = floor(m1/2);
Y = YT(1+nl:n+nl,1+ml:m+ml);

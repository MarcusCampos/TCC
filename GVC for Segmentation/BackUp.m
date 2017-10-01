
%--read the input image
[filename,pathname] = uigetfile('*.jpg','open file');
fname = fullfile(pathname,filename);
%file_name = table(f_name);
file = dir(fullfile(pathname,'*.jpg'));
NF = length(file);
images = cell(NF,1);
for k = 1 : NF
  images{k} = imread(fullfile(pathname, file(k).name));

G = images{k};

[row,col,t] = size(G);
if t>1
G = G(:,:,1);
end

%---preprocessing (filtering,dilating,erosion,etc)--------
G = double(G);             
I = 1 - G/255;          
Ima = Gaus_filter(I,.5);
 
%-----------edge map------------
[Ima_x, Ima_y] = gradient(Ima);
f = sqrt(Ima_x.^2 + Ima_y.^2);
%figure(p); 
imdisp(G); title('original image'); 
     
%------- Compute the external force filed  
[u,v] = MyGVC(f, 30,30,.8,4.6); 
%Nomalizing the GVF external force
mag = sqrt(u.*u+v.*v);
px = u./(mag+1e-10); py = v./(mag+1e-10); % avoid devide by 0 error

%% endocardium segmentation

% initial the contour mannually -> Aqui dimensões

%[a, b] = getpts(G);
%r=7;
%a=90;b=84;r=7;
%a=200;b=50;r=7;
%t = 0:0.05:6.28;
%x = a + r*cos(t);
%y = b + r *sin(t);  
     
 [x1,y1]=ginput(1);
 [x2,y2]=ginput(2);
 r=sqrt((x1-x2).^2+(y1-y2).^2);
  t = 0:0.05:6.28;
  x = x1+ r*cos(t);
 y = y1 + r*sin(t);    
[x,y] = snakeinterp(x,y,1,0.2); 
hold on, 
plot([x;x(1,1)],[y;y(1,1)],'r--','LineWidth',2);
hold off,

% snake deformation with circle-shape constraint
for i=1:8        
    [x,y] = snakedeform_Endocardium(x,y,.5,.5,0.8,1,2,px,py,5);
    [x,y] = snakeinterp(x,y,1,0.2);  %
    hold on
    x = x(:); y = y(:);
    plot([x;x(1,1)],[y;y(1,1)],'r','LineWidth',1);
    title(['Deformation in progress,  iter = ' num2str(i*5)])
    pause(0.5);
end
hold on
    x = x(:); y = y(:);
    plot([x;x(1,1)],[y;y(1,1)],'w','LineWidth',2);
hold off
  
%     figure(p+1); 
     imdisp(G);hold on;
     hold on, 
     plot([x;x(1,1)],[y;y(1,1)],'w','LineWidth',2);
     hold off
     title(['Final result,  iter = ' num2str(i*5)]);
     
%%---------------------------------------------
     figure(3); % display the final segmentation results
     %imdisp(G);hold on;
     imdisp(G);
     hold on, 
     plot([x;x(1,1)],[y;y(1,1)],'w','LineWidth',2);
     hold off
   
     title(['Final result,  iter = ' num2str(i*5)]);
   
%% epicardium segmentation
% compute the new edge map for epicardium
f = NewEdgeMap(Ima, x,y); 
% compute the external force filed
[u,v] = MyGVC(f, 30,30,.8,4.6); %
%Nomalizing the GVF external force
mag = sqrt(u.*u+v.*v);
px = u./(mag+1e-10); py = v./(mag+1e-10); % avoid devide by 0 error

X0=x;
Y0=y;
X00=x;
Y00=y;

n = length(X0);
xc = mean(X0);
yc = mean(Y0);

r = sqrt( (X0 - repmat(xc,n,1)).^2.0 + (Y0 - repmat(yc,n,1)).^2.0 );
rbar = mean(r);

hold on
figure(1)
for i=1:10         %maxium iteration 125
    [x,y] = snakedeform_Epicardium(x,y,0.5,0.5,0.4,0.65,1,px,py,X0,Y0,5,rbar); 
    [x,y] = snakeinterp(x,y,2,0.2);  %dmax=2 dmin=0.5

    len=length(x);
    [X0,Y0] = alignEndoToEpi(len,X00,Y00);
    [x,y] = reAlign(x,y);

    snakedisp(x,y,'r') ;
    title(['Deformation in progress,  iter = ' num2str(i*5)]);
    pause(0.5);
end
hold on, 
plot([x;x(1,1)],[y;y(1,1)],'w','LineWidth',2);
hold off 
title(['Final result,  iter = ' num2str(i*5)]);
%-----------------------------------------------
figure(3); 
imdisp(G);
hold on, 
plot([x;x(1,1)],[y;y(1,1)],'w','LineWidth',2);
hold off 
title(['Final result,  iter = ' num2str(i*5)]);
end


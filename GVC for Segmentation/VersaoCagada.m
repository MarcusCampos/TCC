
%--read the input image
[filename,pathname] = uigetfile('*.bmp','open file');
fname = fullfile(pathname,filename);
%file_name = table(f_name);
file = dir(fullfile(pathname,'*.bmp'));
NF = length(file);
images = cell(NF,1);
for k = 1 : 2
  images{k} = imread(fullfile(pathname, file(k).name));

G = images{k};
%G = imcrop(G, [714.5 0.5 452 456]);

%imdisp(G);
[row,col,t] = size(G);
if t>1
G = G(:,:,1);
end
 imdisp(G);
%---preprocessing (filtering,dilating,erosion,etc)--------
G = double(G);             
I = 1 - G/255;          
Ima = Gaus_filter(I,.5);
 imdisp(Ima);
%-----------edge map------------
[Ima_x, Ima_y] = gradient(Ima);
f = sqrt(Ima_x.^2 + Ima_y.^2);
imdisp(G); title('original image'); 
     
%------- Compute the external force filed  
[u,v] = MyGVC(f, 30,30,.8,4.6); 
%Nomalizing the GVF external force
mag = sqrt(u.*u+v.*v);
px = u./(mag+1e-10); py = v./(mag+1e-10); % avoid devide by 0 error

global XXa ;
global YYa ;

 % alterei aqui o de .2 para 0.8
 if k == 1
    [x1,y1]=ginput(1);
    [x2,y2]=ginput(1);
    t = 0:0.05:6.28; 
    r=sqrt((x1-x2).^2+(y1-y2).^2);
    x = x1+ r*cos(t);
    y = y1 + r*sin(t);
    XXa =x;
    TTa = t;
    YYa=y;
    [x,y] = snakeinterp(x,y,1,0.2);
 else
      %[x,y] = snakeinterp(XXa,YYa,1,0.8);
      [x,y] = snakeinterp(XXa,YYa,1.5,0.05);
end

hold on, 
plot([x;x(1,1)],[y;y(1,1)],'r--','LineWidth',2);
hold off,

% snake deformation with circle-shape constraint
%for i=1:8        
 %   [x,y] = snakedeform_Endocardium(x,y,.5,.5,1,1,2,px,py,5);
for i=1:12        
    [x,y] = snakedeform_Endocardium(x,y,.5,.8,0.8,1,2,px,py,5);

    [x,y] = snakeinterp(x,y,1,0.2);  %
    hold on
    x = x(:); y = y(:);
    plot([x;x(1,1)],[y;y(1,1)],'r','LineWidth',1);
    title(['Deformation in progress,  iter = ' num2str(i*5)])

end
hold off;
hold on
    x = x(:); y = y(:);
    plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',2, 'MarkerEdgeColor','b');
hold off

     figure(k); 
     imdisp(G);
     hold on, 
     plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',2, 'MarkerEdgeColor','b');
     title(['Endocardio Final Segmentation Result,  iter = ' num2str(i*5)]);
     fig = figure(k); 
     hold off;
   
%% epicardium segmentation
% compute the new edge map for epicardium
f = NewEdgeMap(Ima, x,y); 
% compute the external force filed
[u,v] = MyGVC(f, 30,30,.8,4); %
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

figure(1)
tempx = x; %variavel para printar imagens finais - Endocardio
tempy = y; %variavel para printar imagens finais - Endocardio
for i=1:15        %maxium iteration 125
    %snakedeform_Endocardium(x,y,.5,.5,0.8,1,2,px,py,5);
    [x,y] = snakedeform_Epicardium(x,y,0.5,0.3,1,0.1,1,px,py,X0,Y0,5,rbar);
    [x,y] = snakeinterp(x,y,1,0.2);  %dmax=2 dmin=0.5

    len=length(x+y);
    [X0,Y0] = alignEndoToEpi(len,X00,Y00);
    [x,y] = reAlign(x,y);
%for i=1:10         %maxium iteration 125
%    [x,y] = snakedeform_Epicardium(x,y,0.5,0.5,0.4,0.65,1,px,py,X0,Y0,5,rbar); 
%    [x,y] = snakeinterp(x,y,2,0.2);  %dmax=2 dmin=0.5

%    len=length(x);
 %   [X0,Y0] = alignEndoToEpi(len,X00,Y00);
 %   [x,y] = reAlign(x,y);
    snakedisp(x,y,'r') ;
    title(['Deformation in progress,  iter = ' num2str(i*5)]);
    pause(0.5);
end
hold off;
hold on, 
plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',2);
hold off 
title(['Final result,  iter = ' num2str(i*5)]);
    figure(k); 
     imdisp(G);
     hold on, 
     plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',2, 'MarkerEdgeColor','b');
     title(['Enpicardio Final Segmentation Result,  iter = ' num2str(i*5)]);
     fig = figure(k); 
     hold off;
%-----------------------------------------------

figure(k); 
imdisp(G);
hold on, 
title(['Endocardio Final Segmentation Result,  iter = ' num2str(i*5)]);
plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
strgeral = 'EndoCardioSerg';
stringnumer = num2str(k);
fileNameFormat = '.jpg';
stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
saveas(fig, stringNameFile);
hold off 


figure(k); 
imdisp(G);
hold on, 
plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
title(['Enpicardio Final Segmentation Result,  iter = ' num2str(i*5)]);
strgeral = 'EpiCardioSerg';
stringnumer = num2str(k);
fileNameFormat = '.jpg';
stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
saveas(fig, stringNameFile);
hold off 

figure(k); 
imdisp(G);
hold on, 
plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
title(['Final Segmentation Result,  iter = ' num2str(i*5)]);
strgeral = 'Seg';
stringnumer = num2str(k);
fileNameFormat = '.jpg';
stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
saveas(fig, stringNameFile);
hold off 
end


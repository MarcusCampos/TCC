
%--read the input image
[filename,pathname] = uigetfile('*.jpg','open file');
fname = fullfile(pathname,filename);
%file_name = table(f_name);
file = dir(fullfile(pathname,'*.jpg'));
NF = length(file);
images = cell(NF,1);

for k = 1 : 2
    filechoose = randi(NF,1);
  images{k} = imread(fullfile(pathname, file(k).name));

G = images{k};
%G = imcrop(G, [714.5 0.5 452 456]);

%imdisp(G);
[row,col,t] = size(G);
if t>1
G = G(:,:,1);
end
%figure(G);
%---preprocessing (filtering,dilating,erosion,etc)--------
G = double(G);             
I = 1 - G/255;          
Ima = Gaus_filter(I,.5);
if (k== 1)Ima = imgaussfilt(I, 15);end
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

global centerX;
global centerY;
global centerXFirst;
global centerYFirst;
global endoxFirst;
global endoyFirst;
global epixFirst;
global epiyFirst;

 % alterei aqui o de .2 para 0.8
 %if k == 1
    [x1,y1]=ginput(1);
    [x2,y2]=ginput(1);
    [x3,y3]=ginput(1);
    centerX = x1;
    centerY = y1;
    t = 0:0.05:6.28; 
    r=sqrt((x1-x2).^2+(y1-y2).^2);
    x = x1+ r*cos(t);
    y = y1 + r*sin(t);
    XXa =x;
    TTa = t;
    YYa=y;
    [x,y] = snakeinterp(x,y,1,0.2);
%  else
%       [x,y] = snakeinterp(XXa,YYa,1.5,0.05);
% end

hold on, 
plot([x;x(1,1)],[y;y(1,1)],'r--','LineWidth',2);
hold off,

% snake deformation with circle-shape constraint
for i=1:12        
    [x,y] = snakedeform_Endocardium(x,y,.5,.8,0.8,1,2,px,py,5);
    [x,y] = snakeinterp(x,y,1,0.2);  %
    hold on
    x = x(:); y = y(:);
    plot([x;x(1,1)],[y;y(1,1)],'r','LineWidth',1);
    title(['Deformation in progress,  iter = ' num2str(i*5)])

end
hold off;


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

figure(1)
tempx = x; %variavel para printar imagens finais - Endocardio
tempy = y; %variavel para printar imagens finais - Endocardio
control = 0;
for i=1:10 
    %maxium iteration 125
         if i>1 
             dify = tempy - y;
             difx = tempx - x; 
              if all(dify < 2)
                  control = control +1;
              end
         end
             if control >2
                t = 0:0.05:6.28; 
                r=sqrt((x1-x3).^2+(y1-y3).^2);
                x = x1+ r*cos(t);
                y = y1 + r*sin(t);
                [x,y] = snakeinterp(x,y,1,0.2);
                hold on, 
                plot([x;x(1,1)],[y;y(1,1)],'r--','LineWidth',2);
                hold off,

                % snake deformation with circle-shape constraint
                for i=1:12        
                    [x,y] = snakedeform_Endocardium(x,y,.5,.8,0.8,1,2,px,py,5);
                    [x,y] = snakeinterp(x,y,1,0.2);  %
                    hold on
                    x = x(:); y = y(:);
                    plot([x;x(1,1)],[y;y(1,1)],'r','LineWidth',1);
                    title(['Deformation in progress,  iter = ' num2str(i*5)])

                end
                 hold off;
                break;
             else 
                [x,y] = snakedeform_Epicardium(x,y,0.5,0.5,0.4,0.65,1,px,py,X0,Y0,5,rbar); 
                [x,y] = snakeinterp(x,y,2,0.2);  %dmax=2 dmin=0.5
                len=length(x);
                [X0,Y0] = alignEndoToEpi(len,X00,Y00);
                [x,y] = reAlign(x,y);

                snakedisp(x,y,'r') ;
                hold on
                x = x(:); y = y(:);
                plot([x;x(1,1)],[y;y(1,1)],'r','LineWidth',1);
                title(['Deformation in progress,  iter = ' num2str(i*5)]);
            end
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
%saveas(fig, fullfile(pathname, stringNameFile));
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
%saveas(fig,fullfile(pathname, stringNameFile));
%for i = 1: sizeXendo
 %  plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
%    diffX = tempx(i,1) - centerX ;
   % dX = centerX:+.0001:tempx(i,1);
 %   A = [centerX,tempx(i,1)];
    %[diffeX,diffeY] = [A,B];
   %b = linspace(A,B, diffX);
%   dX = A:B;
%   a = linspace(centerY,tempy(i,1), diffY);
%   plot(A,B,'g','LineWidth',1);
%end
%title(['Final Segmentation Result - Selection Endocardio Preenchido,  iter = ' num2str(i*5)]);
%strgeral = 'SegSelection';
%stringnumer = /num2str(k);
%fileNameFormat = '.jpg';
%stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
%saveas(fig, fullfile(pathname, stringNameFile));
%hold off;

figure(k); 
imdisp(G);
hold on, 
plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
% plot(x,y,'y','LineWidth',1);
% rad=sqrt((x1-x(1,1)).^2+(y1-y(1,1)).^2);
% viscircles([x1,y1],rad);
plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
title(['Final Segmentation Result,  iter = ' num2str(i*5)]);
strgeral = 'Seg';
stringnumer = num2str(k);
fileNameFormat = '.jpg';
stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
%saveas(fig, fullfile(pathname, stringNameFile));

% [filename,pathname] = uigetfile('*.jpg','open file');
% fname = fullfile(pathname,filename);
% file = dir(fullfile(pathname,'*.jpg'));
% A = imread(fullfile(pathname, filename));
% imcrop(A, [714.5 0.5 452 456]);
% [filename,pathname] = uigetfile('*.jpg','open file');
% fname = fullfile(pathname,filename);
% file = dir(fullfile(pathname,'*.jpg'));
% B = imread(fullfile(pathname, filename));
% C = imfuse(A,B,'blend','Scaling','joint');
% imshow(C);
hold off ;

 if k == 1
     %
%     limitx = max(x);
%     limiy = max(y);
%     minx= min(x);
%     miny = min(y);
%     diffx = limitx - x;
%     diffy = limiy - y;
%     newvectorX = minx:+1:limitx;
%     newvectorY = miny:+1:limiy;
%     plot([newvectorX;newvectorX(1,:)],[newvectorY;newvectorY(1,:)],'g','LineWidth',1);
%     for i= 0:length(x)
%         NumberNewPoints = 3;
%         xvals = linspace(x(i,1), x(i+1,1), NumberNewPoints+2);
%         yvals = linspace(y(i,1), y(i,1), NumberNewPoints+2);
%         pts = [xvals(:), yvals(:)];
%     end
%     point1=[x,y];
%     point2=[5 5 5];
%     t=0:.01:1;
%     C=repmat(point1,length(t),1)'+(point2-point1)'*t
%     plot([x;x(1,1)],[y;y(1,1)],'g','LineWidth',1);
%     hold off;
%     figureNewp = figure; 
%     hold on,
%     index = miny
%     newvectorX = minx:+1:limitx;
%     newvectorY = miny:+1:limiy;
%     %plot([newvectorX;newvectorX(1,:)],[newvectorY;newvectorY(1,:)],'g','LineWidth',1);
     endoxFirst = tempx;
     endoyFirst = tempy; 
     epixFirst = x; 
     epiyFirst =y;
     centerXFirst = centerX;
     centerYFirst = centerY;
     %hold off;
%  else 
%       CompatorOfSegments(endoxFirst, endoyFirst,epixFirst, epiyFirst,tempx, tempy, x,y)
 end 
 if k == 2
    figureNew = figure; 
    hold on, 
    sizeXendo= size(tempx);
    limitx = max(x);
    limiy = max(y);
    minx= min(x);
    miny = min(y);
    %plot([x;x(1,1)],[y;y(1,1)],'g','LineWidth',1);
  %  plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
    plot([endoxFirst;endoxFirst(1,1)],[endoyFirst;endoyFirst(1,1)],'b','LineWidth',1);
    %plot([epixFirst;epixFirst(1,1)],[epiyFirst;epiyFirst(1,1)],'b','LineWidth',1);
    title(['Final Segmentation Result - Selection Compare,  iter = ' num2str(i*5)]);
    strgeral = 'Salvame';
    stringnumer = num2str(k);
    fileNameFormat = '.jpg';
    stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
    saveas(figureNew, fullfile(pathname, stringNameFile));
    for i =1:size(endoxFirst)
        for k =size(endoyFirst):-1:1 

               A = endoxFirst(i,1) ;
               B = endoyFirst(i,1);
               Al = endoxFirst(k,1) ;
               Bl = endoyFirst(k,1);
               
               for ik =A:Al
                for ki =B : Bl
                   scatter((A),(B));
                   A = A+1
                   B = B+1
                end
              end 
            
        end
    end
    N = verifyIntersection(endoxFirst, endoyFirst, centerXFirst, centerYFirst, x, y, centerX, centerY);
    hold off;
 end
end

A_clr = 'D:\pasa emergencial\Documentos\Faculdade\6°Semestre\7° Semestre\MICCAI\Algoritmo Segment\YUWEI algoritmo\Comparetor\SegSelection1.jpg'; % random image
    B_clr = 'D:\pasa emergencial\Documentos\Faculdade\6°Semestre\7° Semestre\MICCAI\Algoritmo Segment\YUWEI algoritmo\Comparetor\SegSelection1.jpg'; % random image
    A_bw = A_clr>thresh; % convert to B&W
    B_bw = B_clr>thresh; % convert to B&W
    idx = A_bw==B_bw;    % compare
    idy = all(idx,3);    % R==G==B
    out = sum(idy(:));   % count
% figure1 = figure;
% ax1 = axes('Parent',figure1);
% ax2 = axes('Parent',figure1);
% set(ax1,'Visible','off');
% set(ax2,'Visible','off');
% [a,map,alpha] = imread('D:\pasa emergencial\Documentos\Faculdade\6°Semestre\7° Semestre\MICCAI\Algoritmo Segment\YUWEI algoritmo\Comparetor\SegSelection1.jpg');
% I = imshow(a,'Parent',ax2);
% set(I,'AlphaData',alpha);
% imshow('D:\pasa emergencial\Documentos\Faculdade\6°Semestre\7° Semestre\MICCAI\Algoritmo Segment\YUWEI algoritmo\Comparetor\Seg2.jpg','Parent',ax1);
%{
         if i>1 & control==1
         dify = tempy - y;
             difx = tempx - x; 
             if all(dify < 2)
                 [x1,y1]=ginput(1);
                 [x2,y2]=ginput(1);
                 t = 0:0.05:6.28;
                 r=sqrt((x1-x2).^2+(y1-y2).^2);
                 newSelectionX = x1+ r*cos(t);
                 newSelectionY = y1 + r*sin(t);
                 subx = newSelectionX;
                 suby = newSelectionY;
                 x = x(1:340, 1);
                 y = y(1:340, 1);

                [x,y] = snakeinterp(x,y,1,0.2);

                f = NewEdgeMap(Ima, x,y); 
                [u,v] = MyGVC(f, 30,30,.8,4.6); 
                mag = sqrt(u.*u+v.*v);
                px = u./(mag+1e-10); py = v./(mag+1e-10); =
                X0=x;
                Y0=y;

                n = length(X0);
                xc = mean(X0);
                yc = mean(Y0);

                r = sqrt( (X0 - repmat(xc,n,1)).^2.0 + (Y0 - repmat(yc,n,1)).^2.0 );
                rbar = mean(r);
             end
        end
    %}


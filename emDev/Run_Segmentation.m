
%--read the input image

path = 'C:\Users\marcu\Documents\Faculdade\6°Semestre\7° Semestre\MICCAI\Algoritmo Segment\YUWEI algoritmo\Automação';
DIRs = dir(path);
ListagemdeDiretorios = struct2cell(DIRs);

%[filename,pathname] = uigetfile('*.jpg','open file');
%fname = fullfile(pathname,filename);
%file_name = table(f_name);
%file = dir(fullfile(pathname,'*.jpg'));
file = dir(strcat(path,'\',ListagemdeDiretorios{1,3}));
ListagemdeDiretoriosSEMSEG = struct2cell(file);
DiretorioAtual = dir(strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,3}));
ListagemdeDiretorioAtual = struct2cell(DiretorioAtual);
NF = length(DiretorioAtual);
images = cell(NF,1);

for k = 1 : NF
    
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

    if k ==1 
        images{k} = imread((strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,3}, '\',ListagemdeDiretorioAtual{1,k+2})));
        %funcao para setar os pontos de seleção 
        [P1,P2,P3] = getPontosIniciais((strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,3}, '\',ListagemdeDiretorioAtual{1,k+2})));
    else
        images{k} = imread((strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,3}, '\',ListagemdeDiretorioAtual{1,k+2})));
    end


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
%Ima = imgaussfilt(I, 15);
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

 % alterei aqui o de .2 para 0.8
 if k == 1
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
  else
       [x,y] = snakeinterp(XXa,YYa,1.5,0.05);
 end

hold on, 
plot([x;x(1,1)],[y;y(1,1)],'r--','LineWidth',2);
hold off,

% snake deformation with circle-shape constraint
for i=1:20        
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
for i=1:120 
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
saveas(fig, fullfile(pathname, stringNameFile));
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
saveas(fig,fullfile(pathname, stringNameFile));


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
saveas(fig, fullfile(pathname, stringNameFile));
hold off ;

%  if k == 1
%      endoxFirst = tempx;
%      endoyFirst = tempy; 
%      epixFirst = x; 
%      epiyFirst =y;
%      centerXFirst = centerX;
%      centerYFirst = centerY;
%  end 
%  if k == 2
%     figureNew = figure; 
%     hold on, 
%     %imshow(figureNew);
%     %plot([x;x(1,1)],[y;y(1,1)],'g','LineWidth',1);
%   %  plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
%     plot([endoxFirst;endoxFirst(1,1)],[endoyFirst;endoyFirst(1,1)],'k','LineWidth',1);
%     %plot([epixFirst;epixFirst(1,1)],[epiyFirst;epiyFirst(1,1)],'b','LineWidth',1);
%     %title(['Final Segmentation Result - Selection Compare,  iter = ' num2str(i*5)]);
%     axis off
%     hold off;
%     strgeral = 'Salvame';
%     stringnumer = num2str(k);
%     fileNameFormat = '.jpg';
%     stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
%     fullname = fullfile(pathname, stringNameFile);
%     saveas(figureNew, fullfile(pathname, stringNameFile));
%     figureNew2 = figure; 
% 
%     hold on, 
%     plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'k','LineWidth',1);
%      axis off 
%     hold off;
%     strgeral = 'Salvame2';
%     stringnumer = num2str(k);
%     fileNameFormat = '.jpg';
%     stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
%     fullname2 = fullfile(pathname, stringNameFile);
%     saveas(figureNew2, fullfile(pathname, stringNameFile));
%      N = verifyIntersection(fullname,fullname2);
%     hold off;
%  end
end



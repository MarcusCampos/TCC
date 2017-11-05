
%--read the input image

path = 'C:\Users\marcu\Documents\Faculdade\6°Semestre\7° Semestre\MICCAI\Algoritmo Segment\YUWEI algoritmo\Automação';
DIRs = dir(path);
ListagemdeDiretorios = struct2cell(DIRs);

%[filename,pathname] = uigetfile('*.jpg','open file');
%fname = fullfile(pathname,filename);
%file_name = table(f_name);
%file = dir(fullfile(pathname,'*.jpg'));
arquivoSemSegmentacao = dir(strcat(path,'\',ListagemdeDiretorios{1,3}));
arquivoSegmentacao = dir(strcat(path,'\',ListagemdeDiretorios{1,4}));
ListagemdeDiretoriosSEMSEG = struct2cell(arquivoSemSegmentacao);
ListagemdeDiretoriosComSEG = struct2cell(arquivoSegmentacao);
ResultadosSeg = strcat(path,'\Resultados');
Seed = strcat(path,'\Seed');
OkFolder = mkdir(ResultadosSeg)
OkFolder = mkdir(Seed)
if(size(ListagemdeDiretoriosSEMSEG) == size(ListagemdeDiretoriosComSEG))
    disp('pasta mesmo tamanho ')
end
for fileSize = 3:size(ListagemdeDiretoriosSEMSEG)
    DiretorioAtualSemSeg = dir(strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,fileSize}));
    DiretorioAtualComSeg = dir(strcat(path,'\',ListagemdeDiretorios{1,4},'\',ListagemdeDiretoriosComSEG{1,fileSize}));
    
    
    SeedImage = strcat(Seed,'\', ListagemdeDiretoriosComSEG{1,fileSize});
    OkFolderResultadosNSeg = mkdir(SeedImage);
    for i= 3:size(DiretorioAtualSemSeg)
               ArrayBMP = strsplit(DiretorioAtualSemSeg(i,1).name, '.');
               if strcmp(ArrayBMP{2},'bmp') ==1 
                   ArquivoBMP = strcat(ArrayBMP{1}, '.', ArrayBMP{2});
               end 

    end
    movefile((strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,fileSize}, '\',ArquivoBMP)),SeedImage);
    DiretorioAtualSemSeg = dir(strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,fileSize}));
    ListagemdeSemSeg = struct2cell(DiretorioAtualSemSeg);
    ListagemdeComSeg = struct2cell(DiretorioAtualComSeg);
    %Trata Ordenação
    % Convert to a matrix
    ListagemdeSemSegfields = fieldnames(DiretorioAtualSemSeg);
    ListagemdeComSegfields = fieldnames(DiretorioAtualComSeg);
    sz = size(ListagemdeSemSeg)  
    sz1 = size(ListagemdeComSeg)  
    ListagemdeSemSeg = reshape(ListagemdeSemSeg, sz(1), []);      % Px(MxN)
    ListagemdeComSeg = reshape(ListagemdeComSeg, sz1(1), []);

    ListagemdeSemSeg = ListagemdeSemSeg';
    ListagemdeComSeg = ListagemdeComSeg'; % (MxN)xP
    ListagemdeSemSeg = sortrows(ListagemdeSemSeg, 1);
    ListagemdeComSeg = sortrows(ListagemdeComSeg, 1);
    
    ListagemdeSemSeg = reshape(ListagemdeSemSeg', sz);
    ListagemdeComSeg = reshape(ListagemdeComSeg', sz1);
    ListagemdeSemSeg = cell2struct(ListagemdeSemSeg, ListagemdeSemSegfields, 1);
    ListagemdeComSeg = cell2struct(ListagemdeComSeg, ListagemdeComSegfields, 1);

    
    ResultadosNSegString = strcat(ResultadosSeg,'\',ListagemdeDiretoriosSEMSEG{1,fileSize},'\',ListagemdeDiretorios{1,3});
    ResultadosSegString = strcat(ResultadosSeg,'\', ListagemdeDiretoriosComSEG{1,fileSize},'\',ListagemdeDiretorios{1,4});
    
    OkFolderResultadosSeg = mkdir(ResultadosNSegString)
    OkFolderResultadosNSeg = mkdir(ResultadosSegString)
    
    NF = length(ListagemdeSemSeg)-2;
    images = cell(NF,1);
    [~,index] = sortrows({ListagemdeComSeg.name}.'); ListagemdeComSeg = ListagemdeComSeg(index(end:-1:1)); clear index
    [~,index] = sortrows({ListagemdeComSeg.name}.'); ListagemdeComSeg = ListagemdeComSeg(index); clear index
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
    k =1;
     if k ==1 
            
            images{k} = imread((strcat(SeedImage, '\',ArquivoBMP)));
            %funcao para setar os pontos de seleção 
            [P1,P2,P3] = getPontosIniciais((strcat(SeedImage, '\',ArquivoBMP)));
           
     end

    for k = 2 : NF
    
    if (k >= size(ListagemdeComSeg) == 1)
        break;
    end;
    
    images{k} = imread((strcat(path,'\',ListagemdeDiretorios{1,3},'\',ListagemdeDiretoriosSEMSEG{1,fileSize}, '\',ListagemdeSemSeg(k+1,1).name)));
    G = images{k};
    [row,col,t] = size(G);
    if col > 1800
        G = imcrop(G, [714.5 0.5 452 456]);
    end;
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
   % imdisp(Ima);
    %-----------edge map------------
    [Ima_x, Ima_y] = gradient(Ima);
    f = sqrt(Ima_x.^2 + Ima_y.^2);
    %imdisp(G); title('original image'); 

    %------- Compute the external force filed  
    [u,v] = MyGVC(f, 30,30,.8,4.6); 
    %Nomalizing the GVF external force
    mag = sqrt(u.*u+v.*v);
    px = u./(mag+1e-10); py = v./(mag+1e-10); % avoid devide by 0 error

     % alterei aqui o de .2 para 0.8
     if k == 2
        centerX = P1{1};
        centerY = P1{2};
        t = 0:0.05:6.28; 
        r=sqrt((P1{1}-P2{1}).^2+(P1{2}-P2{2}).^2);
        x = P1{1}+ r*cos(t);
        y = P1{2} + r*sin(t);
        XXa =x;
        TTa = t;
        YYa=y;
        [x,y] = snakeinterp(x,y,1,0.2);
      else
           [x,y] = snakeinterp(XXa,YYa,1.5,0.05);
     end
    
    imdisp(G);
    hold on, 
    plot([x;x(1,1)],[y;y(1,1)],'r--','LineWidth',2);
    hold off;

    % snake deformation with circle-shape constraint
    for i=1:15        
        [x,y] = snakedeform_Endocardium(x,y,.5,.8,0.8,1,2,px,py,5);
        [x,y] = snakeinterp(x,y,1,0.15);  %
        hold on
        x = x(:); y = y(:);
        plot([x;x(1,1)],[y;y(1,1)],'r','LineWidth',1);
        title(['Deformation in progress,  iter = ' num2str(i*5)])

    end
    hold off;

    %Salvando o resultado de endocardio
    figure(k); 
    imdisp(G);
    hold on, 
    title(['Endocardio Final Segmentation Result,  iter = ' num2str(i*5)]);
    plot([x;x(1,1)],[y;y(1,1)],'g','LineWidth',1);
    ResultadoTratado = 'EndoCardioSegmentado';
    NomeArquivoOriginal = ListagemdeSemSeg(k+2,1).name;
    stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
    fullname = fullfile(ResultadosNSegString, stringNameFile);
    fig = figure(k); 
    saveas(fig, fullname);
    hold off 



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

    figure(k)
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
                    r=sqrt((P1{1}-P3{1}).^2+(P1{2}-P3{2}).^2);
                    x = P1{1}+ r*cos(t);
                    y = P1{2} + r*sin(t);
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

%     figure(k); 
%     imdisp(G);
%     hold on, 
%     title(['Endocardio Final Segmentation Result,  iter = ' num2str(i*5)]);
%     plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
%     strgeral = 'EndoCardioSerg';
%     stringnumer = num2str(k);
%     fileNameFormat = '.jpg';
%     stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
%     saveas(fig, fullfile(ResultadosNSegString, stringNameFile));
%     hold off 


    figure(k); 
    imdisp(G);
    hold on, 
    plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
    title(['Enpicardio Final Segmentation Result,  iter = ' num2str(i*5)]);
    ResultadoTratado = 'EpiCardioSegmentado';
    NomeArquivoOriginal = ListagemdeSemSeg(k+2,1).name;
    stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
    fullname = fullfile(ResultadosNSegString, stringNameFile);
    saveas(fig,fullname);
    hold off;


    figure(k); 
    imdisp(G);
    hold on, 
    plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
    plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
    title(['Final Segmentation Result,  iter = ' num2str(i*5)]);
    ResultadoTratado = 'SegmentacaoFinalEpiEndoCardio';
    NomeArquivoOriginal = ListagemdeSemSeg(k+2,1).name;
    stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
    fullname = fullfile(ResultadosNSegString, stringNameFile);
    saveas(fig, fullname);
    hold off ;
    
    endoxFirst = tempx;
    endoyFirst = tempy; 
    epixFirst = x; 
    epiyFirst =y;
    centerXFirst = centerX;
    centerYFirst = centerY;
    
    
    %Aqui Começa a comparação com a segmentada 
    images{k} = imread((strcat(path,'\',ListagemdeDiretorios{1,4},'\',ListagemdeDiretoriosComSEG{1,fileSize}, '\',ListagemdeComSeg(k+1,1).name)));
    G = images{k};
    G = imcrop(G, [714.5 0.5 452 456]);
    clearvars x y 
    %G = imcrop(G, [714.5 0.5 452 456]);
    [row,col,t] = size(G);
    if t>1
    G = G(:,:,1);
    end
    %figure(G);
    %---preprocessing (filtering,dilating,erosion,etc)--------
    G = double(G);             
    I = 1 - G/255;          
    Ima = Gaus_filter(I,.5);
    %imdisp(Ima);
    %-----------edge map------------
    [Ima_x, Ima_y] = gradient(Ima);
    f = sqrt(Ima_x.^2 + Ima_y.^2);
    %imdisp(G); title('original image'); 
%
    %------- Compute the external force filed  
    [u,v] = MyGVC(f, 30,30,.8,4.6); 
    %Nomalizing the GVF external force
    mag = sqrt(u.*u+v.*v);
    px = u./(mag+1e-10); py = v./(mag+1e-10); % avoid devide by 0 error
    
    [x,y] = snakeinterp(XXa,YYa,1.5,0.05);


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
                    r=sqrt((P1{1}-P3{1}).^2+(P1{2}-P3{2}).^2);
                    x = P1{1}+ r*cos(t);
                    y = P1{2} + r*sin(t);
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
%     title(['Final result,  iter = ' num2str(i*5)]);
%         figure(k); 
%          imdisp(G);
%          hold on, 
%          plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',2, 'MarkerEdgeColor','b');
%          title(['Enpicardio Final Segmentation Result,  iter = ' num2str(i*5)]);
%          fig = figure(k); 
%          hold off;
    %-----------------------------------------------

    figure(k); 
    imdisp(G);
    hold on, 
    title(['Endocardio Final Segmentation Result,  iter = ' num2str(i*5)]);
    plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'g','LineWidth',1);
    strgeral = 'EndoCardioSeg';
    stringnumer = num2str(k);
    fileNameFormat = '.jpg';
    stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
    saveas(fig, fullfile(ResultadosSegString, stringNameFile));
    hold off 


    figure(k); 
    imdisp(G);
    hold on, 
    plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
    title(['Epicardio Final Segmentation Result,  iter = ' num2str(i*5)]);
    strgeral = 'EpiCardioSerg';
    stringnumer = num2str(k);
    fileNameFormat = '.jpg';
    stringNameFile = strcat(strgeral,stringnumer,fileNameFormat);
    saveas(fig,fullfile(ResultadosSegString, stringNameFile));
    
     figure(k); 
    imdisp(G);
    hold on, 
    plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'b','LineWidth',1);
    plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
    plot([endoxFirst;endoxFirst(1,1)],[endoyFirst;endoyFirst(1,1)],'g','LineWidth',1);
    plot([epixFirst;epixFirst(1,1)],[epiyFirst;epiyFirst(1,1)],'g','LineWidth',1);
    title(['Imagem Final Segmentaçao Ground Truth e Segmentado,  iter = ' num2str(i*5)]);
    ResultadoTratado = 'ImagemFinalMergeGTComSegmentadoSobreImagem';
    NomeArquivoOriginal = ListagemdeSemSeg(k+2,1).name;
    stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
    fullname = fullfile(ResultadosSegString, stringNameFile);
    saveas(fig,fullname);

        
    F1 = figure; 
    hold on, 
       axis off
       plot([endoxFirst;endoxFirst(1,1)],[endoyFirst;endoyFirst(1,1)],'k','LineWidth',1);
       ResultadoTratado = 'EndocardioSoCotorno';
       NomeArquivoOriginal = ListagemdeSemSeg(k+2,1).name;
       stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
       fullname = fullfile(ResultadosSegString, stringNameFile);
       saveas(F1, fullfile(ResultadosSegString, stringNameFile));
    hold off ;
        
    F2 = figure; 
    hold on, 
        plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'k','LineWidth',1);
        axis off
        ResultadoTratado = 'GTEndocardioSoCotorno';
        NomeArquivoOriginal = ListagemdeComSeg(k+1,1).name;
        stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
        fullname2 = fullfile(ResultadosSegString, stringNameFile);
        saveas(F2, fullfile(ResultadosSegString, stringNameFile));
    hold off ;
    
     F3 = figure; 
    hold on, 
        plot([epixFirst;epixFirst(1,1)],[epiyFirst;epiyFirst(1,1)],'k','LineWidth',1);
        axis off
        ResultadoTratado = 'EpicardioSoCotorno';
        NomeArquivoOriginal = ListagemdeComSeg(k+1,1).name;
        stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
        fullname3 = fullfile(ResultadosSegString, stringNameFile);
        saveas(F3, fullfile(ResultadosSegString, stringNameFile));
    hold off ;
    
     F4 = figure; 
    hold on, 
        plot([x;x(1,1)],[y;y(1,1)],'k','LineWidth',1);
        axis off
        ResultadoTratado = 'GTEpicardioSoCotorno';
        NomeArquivoOriginal = ListagemdeComSeg(k+1,1).name;
        stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
        fullname4 = fullfile(ResultadosSegString, stringNameFile);
        saveas(F4, fullfile(ResultadosSegString, stringNameFile));
    hold off ;
    
     F5 = figure; 
    hold on, 
        plot([tempx;tempx(1,1)],[tempy;tempy(1,1)],'b','LineWidth',1);
        plot([x;x(1,1)],[y;y(1,1)],'b','LineWidth',1);
        plot([endoxFirst;endoxFirst(1,1)],[endoyFirst;endoyFirst(1,1)],'g','LineWidth',1);
        plot([epixFirst;epixFirst(1,1)],[epiyFirst;epiyFirst(1,1)],'g','LineWidth',1);
        axis off
        ResultadoTratado = 'ImagemFinalMergeGTComSegmentadoContorno';
        NomeArquivoOriginal = ListagemdeComSeg(k+1,1).name;
        stringNameFile = strcat(ResultadoTratado,NomeArquivoOriginal);
        fullname5 = fullfile(ResultadosSegString, stringNameFile);
        saveas(F5, fullfile(ResultadosSegString, stringNameFile));
    hold off ;
        
        
          N = verifyIntersection(fullname,fullname2);
          N2 = verifyIntersection(fullname3,fullname4);
        hold off;
     

    end

end 



